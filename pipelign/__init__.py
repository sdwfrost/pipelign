"""
A pipeline for automated alignment
"""

import sys, os, shutil, subprocess, argparse, textwrap
import copy
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO, Phylo
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.motifs import Motif
import tempfile
import time
import joblib
import hashlib
from ete4 import Tree

NUCLEOTIDE_CONSENSUS_ALPHABET = "ACGTURYSWKMBDHVN"
PROTEIN_CONSENSUS_ALPHABET = "ACDEFGHIKLMNPQRSTVWYBXZJUO"
_RESUME_EXISTING = False
_RESUME_CHECKSUMS = {}


class PipelignError(RuntimeError):
    """Raised when the pipeline cannot safely continue."""


EXTERNAL_TOOL_GROUPS = {
    "common": [("mafft",), ("hmmbuild",), ("hmmpress",), ("iqtree", "iqtree2")],
    "dna_cluster": [("cd-hit-est",)],
    "aa_cluster": [("cd-hit",)],
    "dna_search": [("makeblastdb",), ("blastn",), ("nhmmscan",)],
    "aa_search": [("makeblastdb",), ("blastp",), ("hmmscan",)],
}

EXTERNAL_TOOL_DETAILS = {
    ("mafft",): ("mafft", "sequence alignment and profile merging", "stages 1-5"),
    ("hmmbuild",): ("hmmer", "building cluster HMMs", "stages 1-5"),
    ("hmmpress",): ("hmmer", "indexing the HMM database", "stages 1-5"),
    ("iqtree", "iqtree2"): ("iqtree", "cluster-representative tree inference", "stages 1-5"),
    ("cd-hit-est",): ("cd-hit", "DNA/RNA sequence clustering", "stages 1-5 for dna/rna"),
    ("cd-hit",): ("cd-hit", "protein sequence clustering", "stages 1-5 for aa"),
    ("makeblastdb",): ("blast", "BLAST database construction", "stages 3-5"),
    ("blastn",): ("blast", "DNA/RNA fragment search", "stages 3-5 for dna/rna"),
    ("blastp",): ("blast", "protein fragment search", "stages 3-5 for aa"),
    ("nhmmscan",): ("hmmer", "DNA/RNA HMM fragment search", "stages 3-5 for dna/rna"),
    ("hmmscan",): ("hmmer", "protein HMM fragment search", "stages 3-5 for aa"),
}


def _first_available_executable(*names):
    for name in names:
        resolved = shutil.which(name)
        if resolved:
            return resolved
    return None


def _format_tool_names(names):
    return " or ".join(names)


def external_tool_details(names):
    package, purpose, stages = EXTERNAL_TOOL_DETAILS.get(names, ("unknown", "pipelign execution", "selected stages"))
    return {
        "names": names,
        "display": _format_tool_names(names),
        "package": package,
        "purpose": purpose,
        "stages": stages,
    }


def required_external_tools(alphabet="dna", run="G", stage=5):
    tools = list(EXTERNAL_TOOL_GROUPS["common"])
    tools.extend(EXTERNAL_TOOL_GROUPS["aa_cluster" if alphabet == "aa" else "dna_cluster"])

    if stage in (3, 4, 5):
        tools.extend(EXTERNAL_TOOL_GROUPS["aa_search" if alphabet == "aa" else "dna_search"])

    seen = set()
    unique = []
    for names in tools:
        if names not in seen:
            seen.add(names)
            unique.append(names)
    return unique


def missing_external_tools(alphabet="dna", run="G", stage=5):
    return [names for names in required_external_tools(alphabet, run, stage) if _first_available_executable(*names) is None]


def check_external_dependencies(alphabet="dna", run="G", stage=5):
    missing = missing_external_tools(alphabet, run, stage)
    if not missing:
        return

    msg = "\nERROR: Missing required external executable(s):\n"
    for names in missing:
        details = external_tool_details(names)
        msg += "  - {display} ({purpose}; install package: {package})\n".format(**details)
    msg += "Install the missing tools with conda/micromamba/pixi from bioconda/conda-forge, then rerun pipelign.\n"
    raise PipelignError(msg)


def _parallel_worker_plan(task_count, thread_count):
    if task_count < 1:
        raise ValueError("task_count must be positive")
    if thread_count < 1:
        raise ValueError("thread_count must be positive")

    jobs = min(task_count, thread_count)
    return jobs, max(1, thread_count // jobs)


def _records_by_unique_id(records, source):
    records_by_id = {}
    for record in records:
        if not record.id:
            raise ValueError("empty FASTA identifier in <%s>" % source)
        if record.id in records_by_id:
            raise ValueError("duplicate FASTA identifier %r in <%s>" % (record.id, source))
        records_by_id[record.id] = record
    return records_by_id


def _fields(line, minimum, source, line_number):
    fields = line.split()
    if len(fields) < minimum:
        raise ValueError(
            "malformed %s output at line %d: expected at least %d fields"
            % (source, line_number, minimum)
        )
    return fields


def _cluster_number_from_name(name, source):
    parts = name.split(".")
    if len(parts) < 2 or not parts[1].isdigit():
        raise ValueError("could not parse cluster number from %r in %s output" % (name, source))
    return parts[1]


def _validate_text_artifact(path, label):
    if not os.path.isfile(path):
        raise RuntimeError("%s did not produce <%s>" % (label, path))
    if os.path.getsize(path) == 0:
        raise RuntimeError("%s did not produce a nonempty <%s>" % (label, path))


def _validate_hmm_artifact(path, label):
    _validate_text_artifact(path, label)
    with open(path, encoding="utf-8", errors="replace") as handle:
        first_line = handle.readline().strip()
    if not first_line.startswith("HMMER"):
        raise RuntimeError("%s produced an invalid HMM file <%s>" % (label, path))


def _validate_fasta_artifact(
    path,
    label,
    alignment=False,
    expected_ids=None,
    expected_ungapped=None,
):
    _validate_text_artifact(path, label)
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise RuntimeError("%s produced an unreadable FASTA file <%s>" % (label, path))
    records_by_id = _records_by_unique_id(records, path)
    if expected_ids is not None and set(records_by_id) != set(expected_ids):
        raise RuntimeError(
            "%s output identifiers do not match its input identifiers" % label
        )
    if expected_ungapped is not None:
        observed = {
            record.id: str(record.seq).replace("-", "").replace(".", "").upper()
            for record in records
        }
        if observed != expected_ungapped:
            raise RuntimeError(
                "%s output does not preserve its input residues" % label
            )
    if alignment and len({len(record.seq) for record in records}) != 1:
        raise RuntimeError("%s produced sequences with inconsistent aligned lengths" % label)
    return records


def _record_sequences(records):
    return {
        record.id: str(record.seq).replace("-", "").replace(".", "").upper()
        for record in records
    }


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def configure_resume(enabled, checksums=None):
    """Enable reuse of already validated external-tool artifacts."""

    global _RESUME_EXISTING, _RESUME_CHECKSUMS
    _RESUME_EXISTING = bool(enabled)
    _RESUME_CHECKSUMS = dict(checksums or {})


def _checksum_matches(path):
    if not _RESUME_EXISTING or not os.path.isfile(path):
        return False
    expected = _RESUME_CHECKSUMS.get(os.path.normpath(path))
    return expected is not None and _sha256(path) == expected


def _reusable_file(path):
    return _checksum_matches(path) and os.path.getsize(path) > 0


def _reusable_fasta(
    path,
    alignment=False,
    expected_ids=None,
    expected_ungapped=None,
):
    if not _RESUME_EXISTING:
        return False
    if not _checksum_matches(path):
        return False
    try:
        _validate_fasta_artifact(
            path,
            "resume",
            alignment=alignment,
            expected_ids=expected_ids,
            expected_ungapped=expected_ungapped,
        )
    except (OSError, RuntimeError, ValueError):
        return False
    return True


def _copy_fasta(source, destination, expected_ids=None):
    source_records = _validate_fasta_artifact(source, "copy input")
    shutil.copyfile(source, destination)
    _validate_fasta_artifact(
        destination,
        "copy",
        alignment=source.endswith(".aln"),
        expected_ids=expected_ids,
        expected_ungapped=_record_sequences(source_records),
    )


def _run_threaded(tasks, worker, thread_count):
    if not tasks:
        return
    jobs, threads_per_job = _parallel_worker_plan(len(tasks), thread_count)
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = [executor.submit(worker, task, threads_per_job) for task in tasks]
        for future in futures:
            future.result()


def _node_is_leaf(node):
    is_leaf = getattr(node, "is_leaf")
    return is_leaf() if callable(is_leaf) else bool(is_leaf)


def _read_named_tree(tree_file):
    return Tree(tree_file, parser=1)


def _find_tree_node(tree, name):
    matches = list(tree.search_nodes(name=name))
    if not matches:
        raise ValueError("Tree node %r was not found" % name)
    if len(matches) > 1:
        raise ValueError("Tree node name %r is not unique" % name)
    return matches[0]


def _iter_tree_leaves(tree):
    leaves = getattr(tree, "leaves", None)
    if callable(leaves):
        return list(leaves())
    return [node for node in tree.traverse() if _node_is_leaf(node)]


def consensus_from_alignment(aln, alphabet):
    default_alphabet = PROTEIN_CONSENSUS_ALPHABET if alphabet == "aa" else NUCLEOTIDE_CONSENSUS_ALPHABET
    ambiguous = "X" if alphabet == "aa" else "N"
    records = []
    observed = set(default_alphabet)
    lengths = set()

    for record in aln:
        sequence = str(record.seq).upper()
        lengths.add(len(sequence))
        observed.update(char for char in sequence if char != "-")
        records.append(SeqRecord(Seq(sequence), id=record.id, description=""))

    if not records:
        raise ValueError("cannot calculate consensus from an empty alignment")
    if len(lengths) != 1:
        raise ValueError("cannot calculate consensus from sequences with different lengths")
    aln_len = lengths.pop()
    if aln_len == 0:
        raise ValueError("cannot calculate consensus from an alignment with no columns")
    if all(char == "-" for record in records for char in str(record.seq)):
        return ambiguous * aln_len

    gap_only_columns = {
        idx
        for idx in range(aln_len)
        if all(str(record.seq)[idx] == "-" for record in records)
    }
    if gap_only_columns:
        records = [
            SeqRecord(
                Seq("".join(char for idx, char in enumerate(str(record.seq)) if idx not in gap_only_columns)),
                id=record.id,
                description="",
            )
            for record in records
        ]

    motif_alphabet = "".join(sorted(observed - {"-"}))
    normalized = MultipleSeqAlignment(records)
    consensus = str(Motif(motif_alphabet, normalized.alignment).counts.calculate_consensus(identity=0.5))
    if gap_only_columns:
        consensus_parts = []
        consensus_idx = 0
        for idx in range(aln_len):
            if idx in gap_only_columns:
                consensus_parts.append(ambiguous)
            else:
                consensus_parts.append(consensus[consensus_idx])
                consensus_idx += 1
        consensus = "".join(consensus_parts)
    return consensus.replace("-", ambiguous)

#************************************
class Struct:
    '''
     - this class creates a structure of the parameter values needed for different functions.
     - parameters can be accessed as structName.itemName e.g. mArgs.gCode
    '''
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self,k,v)

#*************************************

class MyStruct(Struct):
    pass

#*****************************************

def lengthThreshold(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError('%r not in range [0.0, 1.0]' % (x,))
    return x

#***************************************
def simThreshold(x):
    x = float(x)
    if x < 0.4 or x > 1.0:
        raise argparse.ArgumentTypeError('%r not in range [0.4, 1.0]' % (x,))
    return x

#***************************************
def iterThreshold(x):
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError('%r must be a non-zero positive integer' % (x,))
    return x

#***************************************
def runType(x):
    x = x.upper()
    if x not in ['J','G']:
        raise argparse.ArgumentTypeError('value of -r must be either J or G')
    return x

#***************************************
def mergeType(x):
    x = x.upper()
    if x not in ['P','C']:
        raise argparse.ArgumentTypeError('value of -e must be either P or C')
    return x

#***************************************

def genCodeLimit(x):
    x = int(x)
    gCodes = [1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25]
    if x not in gCodes:
        raise argparse.ArgumentTypeError('%r is not a valid genetic code' % (x,))
    else:
        return x

#****************************************

def cZip(cDir,tName,zName):
    '''
    creates a zip file of the temporary directory
    '''
    os.chdir(cDir)

    #zName = 'pipelign.' + time.strftime('%Y-%m-%d-%H%M%S')
    try:
        shutil.make_archive(zName,'zip',tName)
    except OSError as e:
        raise PipelignError(str(e)) from e

    print('\nArchive for all temporary files created in %s.zip\n' % zName)
    return

#******************************************

def concatLogFiles(original,newLog):
    '''
    Concatenates two log files
    '''

    lh = open(original,'a')

    with open(newLog,'r') as f:
        for line in f:
            lh.write(line)

    lh.close()

#******************************************
def validateFastaType(fName,alphabet):
    '''
        - Checks the sequence characters in the FASTA file
        - if count(A,T,C,G) > 50%, assumes DNA; assume protein otherwise
        - if alphabet matches with the prediction, returns True
        - otherwise returns False
    '''
    # read in the sequence file and create a string of characters < 2000 length
    residues = ''
    # get around 2000 residues
    try:
        for fasta in SeqIO.parse(fName,'fasta'):
            residues += str(fasta.seq)

            if len(residues) > 2000:
                break
    except FileNotFoundError as e:
        raise PipelignError(str(e)) from e

    # if valid fasta file
    if residues:
        normalized = residues.upper()
        if alphabet in ('dna', 'rna'):
            if alphabet == 'dna' and 'U' in normalized:
                return False
            if alphabet == 'rna' and 'T' in normalized:
                return False
            bases = 'ACGT' if alphabet == 'dna' else 'ACGU'
            nucleotide_count = sum(normalized.count(base) for base in bases)
            return nucleotide_count >= len(residues) * 0.5

        dna_count = sum(normalized.count(base) for base in 'ACGT')
        return dna_count < len(residues) * 0.5
    else:
        raise PipelignError('\nERROR: Could not parse FASTA input file. Please run pipelign with a valid FASTA file\n')

#******************************************

def linsi(seqFile,alnFile,thread,mIterLong,log,p=None):
    '''
      - aligns sequences using MAFFT's L-INS-i method
    '''

    if p:
        lName = log + '.' + str(p)
    else:
        lName = log + '.temp'

    input_records = _validate_fasta_artifact(seqFile, "MAFFT input")
    expected_ids = [record.id for record in input_records]
    expected_ungapped = _record_sequences(input_records)
    if _reusable_fasta(
        alnFile,
        alignment=True,
        expected_ids=expected_ids,
        expected_ungapped=expected_ungapped,
    ):
        return

    cl = ['mafft', '--localpair', '--thread', str(thread), '--maxiterate', str(mIterLong), '--preservecase', seqFile]

    with open(lName, 'a') as lh, open(alnFile, 'w') as ah:
        subprocess.check_call(cl,stdout=ah,stderr=lh)
    _validate_fasta_artifact(
        alnFile,
        "MAFFT",
        alignment=True,
        expected_ids=expected_ids,
        expected_ungapped=expected_ungapped,
    )

    if not p:
        concatLogFiles('pipelign.log',lName)
        os.remove(lName)

#***********************************************************************

def fftnsi(seqFile,alnFile,thread,mIterLong,log,p=None):
    '''
      - aligns sequences using MAFFT's FFT-NS-i method
    '''

    if p:
        lName = log + '.' + str(p)
    else:
        lName = log + '.temp'

    input_records = _validate_fasta_artifact(seqFile, "MAFFT input")
    expected_ids = [record.id for record in input_records]
    expected_ungapped = _record_sequences(input_records)
    if _reusable_fasta(
        alnFile,
        alignment=True,
        expected_ids=expected_ids,
        expected_ungapped=expected_ungapped,
    ):
        return

    cl = ['mafft', '--thread', str(thread), '--maxiterate', str(mIterLong), '--preservecase', seqFile]

    with open(lName, 'a') as lh, open(alnFile, 'w') as ah:
        subprocess.check_call(cl, stdout=ah, stderr=lh)
    _validate_fasta_artifact(
        alnFile,
        "MAFFT",
        alignment=True,
        expected_ids=expected_ids,
        expected_ungapped=expected_ungapped,
    )

    if not p:
        concatLogFiles('pipelign.log',lName)
        os.remove(lName)

#***********************************************************************

def checkPresenceOfFile(fName):
    '''
      - Check whether file <fName> exists;
      - Terminate Pipelign run otherwise
    '''

    if not os.path.exists(fName) or os.stat(fName).st_size == 0:
        return False

    return True

#************************************************************************

def checkMPI():
    cmd = "mpirun"
    fh = open('mpi.log','a')
    mpi = True
    try:
        subprocess.call([cmd],stdout=fh,stderr=fh)
    except argparse.ArgumentTypeError as e:
        print('\nmpirun was not found. Running without MPI\n')
        mpi = False

    fh.close()

    return mpi

#************************************************************************
def makeOutDir(dName,cFlag):
    '''
        - creates output directory to hold pipelign intermediate files
    '''
    if os.path.isdir(dName):
        if cFlag:
            shutil.rmtree(dName)
            print('\nRemoving existing output directory <{}>'.format(dName))
        else:
            msg = '\nERROR: Directory exists with name <{}>. Please use another name for output directory.'.format(dName)
            msg += ' Run pipelign with -c option to overwrite.\n'
            raise PipelignError(msg)

    elif os.path.exists(dName):
        raise PipelignError('\nERROR: File exists with name <{}>. Please use another name for output directory'.format(dName))

    try:
        os.mkdir(dName)
        print('\nCreated output directory <{}>'.format(dName))
    except OSError as e:
        raise PipelignError('\nERROR: system could not create temporary directory. Please try again\n') from e


#**********************************************************************
def makeTempDir(tDir,zName):
    '''
      - creates a temporary directory
      - creates inside user specified directory if path given
    '''

    if tDir is None: # no path provided for temp directory
        try:
            tempDir = tempfile.TemporaryDirectory() # create temporary directory to hold intermediary files
            tDirName = tempDir.name
        except OSError as e:
            raise PipelignError('\nError: system could not create temporary directory. Please try again') from e

    else: # user provided path using -d
        if os.path.exists(tDir):
            tempDir = os.path.join(tDir, zName)
            tDirName = tempDir
            try:
                os.mkdir(tempDir)
            except OSError as e:
                raise PipelignError('\nError: system could not create temporary directory. Please try again') from e
        else:
            raise PipelignError('\nError: Path for temporary directory does not exists. Please run again with correct path.')

    return tempDir, tDirName

#************************************************************************

def copyFile(sName,dName):
    '''
      - This function makes copy of a single file
      - file <sName> is copied to the destination <dName>
    '''

    try:
        shutil.copyfile(sName,dName)
    except OSError as e:
        raise PipelignError(str(e)) from e


def publishAlignment(source, destination):
    """Atomically validate and publish a final alignment."""

    source_records = _validate_fasta_artifact(source, "final alignment", alignment=True)
    destination = os.path.abspath(destination)
    parent = os.path.dirname(destination)
    if not os.path.isdir(parent):
        raise OSError("output directory does not exist: %s" % parent)
    descriptor, temporary = tempfile.mkstemp(
        prefix=".%s-" % os.path.basename(destination),
        suffix=".tmp",
        dir=parent,
    )
    os.close(descriptor)
    try:
        shutil.copyfile(source, temporary)
        _validate_fasta_artifact(
            temporary,
            "final alignment publication",
            alignment=True,
            expected_ids=[record.id for record in source_records],
            expected_ungapped=_record_sequences(source_records),
        )
        with open(temporary, "rb") as handle:
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
        try:
            directory_fd = os.open(parent, os.O_RDONLY)
        except OSError:
            pass
        else:
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)


#************************************************************************

def deAlign(iFile, dFile):
    '''
      - Removes gaps (if any) from the input sequence file
    '''

    #print("\nRemoving the gap characters '-'/'.' from the sequences")

    # Reading sequence file in fasta format
    seqs = list(SeqIO.parse(iFile,'fasta'))

    if len(seqs) >= 1: # at least one sequence present | file is in FASTA format
        _records_by_unique_id(seqs, iFile)

        st = ''

        for seq in seqs:
          st += '>' + seq.id + '\n' + str(seq.seq).replace('-','').replace('.','') + '\n'

        fh = open(dFile,'w')
        fh.write(st)
        fh.close()

        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Gapless sequence file written in <%s>\n' % dFile
        print(msg)

    else: # no sequence present or wrong format
        msg = '\n\nError: Could not read the input sequence file.'
        msg += '\n       Make sure the file is in FASTA format'
        msg += '\n       and at least one sequnce present in file\n'
        raise PipelignError(msg)

#*************************************************************************

def separateFullFragment(seqFile, lenThr, longName, fragName):
    '''
      Reads in the input sequence file.
          - finds length of the longest sequence
          - calculates minimum length required for full length sequences
          - writes full length sequences and fragments into two separate files
    '''

    maxLen = 0

    # get maximum length of the input sequences
    handle = open(seqFile, "r")
    for record in SeqIO.parse(handle, "fasta"):
        if len(record.seq) > maxLen:
            maxLen = len(record.seq)
    handle.close()

    # calculate minimum length for long sequences
    minLengthLong = int(lenThr * maxLen)

    longS = []
    fragS = []

    # create separate lists for long sequences and fragments
    handle = open(seqFile, "r")
    for record in SeqIO.parse(handle, "fasta"):
        if len(record.seq) < minLengthLong:
            fragS.append(record)
        else:
            longS.append(record)
    handle.close()

    # write long file
    SeqIO.write(longS,longName,'fasta')

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' <%s> ' % longName

    fEmpty = True

    if len(fragS) > 0:
        SeqIO.write(fragS, fragName,'fasta')
        msg += ' and <%s> ' % fragName
        fEmpty = False

    msg += 'created\n'
    print(msg)

    return fEmpty, len(fragS)

#************************************************************************
def separateLongAmbigs(alphabet,ambigPer):
    '''
        Separates the long sequences that have ambiguous characters more than a limit
    '''

    # list the ambiguous characters

    if alphabet == 'aa':
        chars = ['X','!','*','?']
    else:
        chars = ['N','X','!','*','?']

    # count ambiguous characters for each sequences
    # if count is more than the threshold, separate them in another file

    gs = [] # contains good quality sequences
    bs = [] # contains bad quality sequences

    handle = open('long.fas','r')

    # counts number of ambigs for each sequences
    for seq in SeqIO.parse(handle,'fasta'):
        seqStr = str(seq.seq).upper()
        acount = 0
        for i in range(len(chars)):
            acount += seqStr.count(chars[i])

        propChars = int(len(seq.seq) * float(ambigPer))

        if acount > propChars:
            bs.append(seq)
        else:
            gs.append(seq)

    if len(gs) > 0:
        SeqIO.write(gs,'long.good.fas','fasta')
    else:
        msg = '\nAll long sequences have ambiguous characters more than allowed proportion %s' % ambigPer
        msg += '\nPlease check whether the correct alphabet (e.g. dna) was selected'
        msg += '\n\nPipelign is exiting\n'
        handle.close()
        raise PipelignError(msg)

    if len(bs) > 0:
        SeqIO.write(bs,'long.bad.fas','fasta')

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Good quality long sequences written in <long.good.fas>\n'
    print(msg)

    handle.close()

#************************************************************************

def runCDHIT(longName,alphabet,per,thread):
    '''
      CD-HIT is used to group similar sequences together in clusters for alignment
    '''

    if _reusable_file('grp') and _reusable_file('grp.clstr'):
        try:
            input_ids = set(
                _records_by_unique_id(
                    list(SeqIO.parse(longName, 'fasta')), longName
                )
            )
            clustered_ids = set()
            with open('grp.clstr') as cluster_handle:
                for line_number, line in enumerate(cluster_handle, 1):
                    if 'Cluster' in line:
                        continue
                    fields = _fields(line, 3, "CD-HIT", line_number)
                    sequence_id = fields[2].replace('>', '').replace('...', '')
                    if sequence_id in clustered_ids:
                        raise ValueError("duplicate CD-HIT assignment for %r" % sequence_id)
                    clustered_ids.add(sequence_id)
            if clustered_ids != input_ids:
                raise ValueError("CD-HIT output does not cover every input sequence")
            return
        except (OSError, RuntimeError, ValueError):
            pass

    # count number of sequences in long sequence file
    seqCount = 0
    handle = open(longName,'r')
    for record in SeqIO.parse(handle,'fasta'):
        seqCount += 1

    '''
    # only one long sequence, one cluster
    if seqCount == 1:
        try:
            shutil.copy(longName,'grp')
        except OSError as e:
            print(e)
            cZip(cDir,tName,zName)
        print('Only one Long sequence present')
        return
    '''
    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'

    # argument string for CD-HIT
    if alphabet == 'dna' or alphabet == 'rna':
        # choose word size for cd-hit-est
        if per >= 0.9:
            ws = 8
        elif per >= 0.88:
            ws = 7
        elif per >= 0.85:
            ws = 6
        elif per >= 0.80:
            ws = 5
        else:
            ws = 4

        cl = ['cd-hit-est','-c', str(per), '-n', str(ws), '-i', longName, '-o', 'grp', '-d', '0', '-T', str(thread)]
        msg += ' CD-HIT-EST started\n'

    elif alphabet == 'aa':
        # choose word size for cd-hit
        if per >= 0.7:
            ws = 5
        elif per >= 0.6:
            ws = 4
        elif per >= 0.5:
            ws = 3
        else:
            ws = 2

        cl = ['cd-hit','-c', str(per), '-n', str(ws), '-i', longName, '-o', 'grp', '-d', '0', '-T', str(thread)]
        msg += ' CD-HIT started\n'

    print(msg)

    with open('cdhit.log', 'a') as lh:
        subprocess.check_call(cl, stdout=lh, stderr=lh)
    _validate_fasta_artifact('grp', "CD-HIT")
    _validate_text_artifact('grp.clstr', "CD-HIT")

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' CD-HIT finished. Files created: <grp> and <grp.clstr>\n'
    print(msg)

    concatLogFiles('pipelign.log','cdhit.log')
    #os.remove('cdhit.log')

#*************************************************************************
def makeClusters(longName,cName):
    '''
      Separate files are created for each clusters
    '''

    if not os.path.exists(cName) or os.stat(cName).st_size == 0:
        msg = '\nError: the file <%s>  could not be found in the given path' % cName
        msg += '\n\nPipelign is exiting\n'
        raise PipelignError(msg)

    # read in the long sequence file
    seqs = list(SeqIO.parse(longName,'fasta'))

    #clsSize = list() # will hold a list of cluster sizes

    # only one long sequence
    if len(seqs) == 1:
        copyFile(longName,'long.0.fas')
        #clsSize.append(1)
        fh = open('long.ClusterList.txt','w')
        fh.write('%s\t0' %seqs[0].id)
        fh.close()
        return 1 #, clsSize

    #cName = 'grp.clstr'

    with open(cName,'r') as fh:
        lines = [line.strip('\n') for line in fh] # read cluster file

    start = 0 # flag for the beginning of first cluster list

    cSeq = [] # hold sequences of a cluster
    assigned_ids = set()

    cls = 0 # count clusters

    st = '' # clusterList string

    seqs_by_id = _records_by_unique_id(seqs, longName)

    # read cluster file and make separate cluster files
    if not lines or 'Cluster' not in lines[0]:
        msg = '\n\nError: <%s> does not contain any cluster' % cName
        msg += '\nPlease try running the program again\n'
        raise PipelignError(msg)

    for i in range(1, len(lines)):
        if 'Cluster' in lines[i]: # start of a new cluster list
            gName = 'long.' + str(cls) + '.fas'
            cls = cls + 1
            if len(cSeq) > 0:
                SeqIO.write(cSeq,gName,'fasta')
                #clsSize.append(len(cSeq))
            cSeq = []

        else: # continue with the existing cluster
            fields = _fields(lines[i], 3, "CD-HIT", i + 1)
            seqID = fields[2].replace('>','').replace('...','')
            if seqID in assigned_ids:
                raise ValueError("duplicate CD-HIT assignment for %r" % seqID)
            assigned_ids.add(seqID)
            st += seqID + '\t' + str(cls) + '\n' # updating clusterList file content
            try:
                cSeq.append(seqs_by_id[seqID])
            except KeyError as exc:
                raise ValueError(
                    "CD-HIT output references unknown sequence identifier %r" % seqID
                ) from exc

    gName = 'long.' + str(cls) + '.fas'
    cls = cls + 1
    if len(cSeq) > 0:
        SeqIO.write(cSeq,gName,'fasta')

    if assigned_ids != set(seqs_by_id):
        missing = sorted(set(seqs_by_id) - assigned_ids)
        raise ValueError(
            "CD-HIT output does not assign every input sequence; missing: %s"
            % ", ".join(missing)
        )

    fh = open('long.ClusterList.txt','w')
    fh.write(st)
    fh.close()

    if cls > 0:
        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' %d cluster file(s) created. File names: <long.x.fas>\n' % cls
        print(msg)
    else:
        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Could not create cluster files.'
        msg += '\nPipelign is exiting\n'
        raise PipelignError(msg)

    return cls #, clsSize

#***********************************************************************
def addClusterNumberToReps(repName,lstFile,outFile):
    '''
        - Reads in the cluster representative FASTA file <grp> and the <long.ClusterList.txt> file
        - Adds cluster number and size to the sequence header e.g. >seq1_size_cluster
        - Temporary file <clsReps.fas> is written
    '''

    # read in the <long.ClusterList.txt> file
    #cList = [line.strip() for line in open(lstFile,'r')]

    clusters_by_id = {}
    cluster_sizes = {}

    with open(lstFile,'r') as f:
        for line_number, line in enumerate(f, 1):
            words = _fields(line, 2, lstFile, line_number)
            if words[0] in clusters_by_id:
                raise ValueError("duplicate sequence identifier %r in <%s>" % (words[0], lstFile))
            clusters_by_id[words[0]] = words[1]
            cluster_sizes[words[1]] = cluster_sizes.get(words[1], 0) + 1

    # read in the cluster reps file
    seqs = list(SeqIO.parse(repName,'fasta'))
    _records_by_unique_id(seqs, repName)

    for seq in seqs:
        if seq.id in clusters_by_id:
            cluster_number = clusters_by_id[seq.id]
            seq.id = seq.id + '_' + str(cluster_sizes[cluster_number]) + '_' + cluster_number
            seq.name = seq.id
            seq.description = seq.id
        else:
            msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
            msg += ' Error: %s was not found in <long.ClusterList.txt>.' % seq.id
            msg += '\nPipelign is exiting\n'
            raise PipelignError(msg)

    # write reps file
    SeqIO.write(seqs,outFile,'fasta')

#***********************************************************************
def makeIQTree(alnFile,thread,alpha):
    '''
      - Constructs phylogenetic tree using IQ-TREE
    '''

    #print('\nCreating IQ-TREE from %s' % alnFile)
    _validate_fasta_artifact(alnFile, "IQ-TREE input", alignment=True)
    tree_path = alnFile + '.treefile'
    if _reusable_file(tree_path):
        try:
            Phylo.read(tree_path, 'newick')
            return
        except Exception:
            pass

    iqtree = _first_available_executable("iqtree", "iqtree2")
    if iqtree is None:
        raise PipelignError("\nERROR: Missing required external executable: iqtree or iqtree2\n")

    if alpha == 'dna' or alpha == 'rna':
        cl = [iqtree, '-s', alnFile, '-m', 'GTR+R4', '-nt', str(thread)]
    elif alpha == 'aa':
        cl = [iqtree, '-s', alnFile, '-m', 'WAG', '-nt', str(thread)]

    with open('tree.log', 'a') as lh:
        subprocess.check_call(cl,stdout=lh,stderr=lh)
    _validate_text_artifact(tree_path, "IQ-TREE")
    try:
        Phylo.read(tree_path, 'newick')
    except Exception as exc:
        raise RuntimeError("IQ-TREE produced invalid Newick output <%s>" % tree_path) from exc

    concatLogFiles('pipelign.log','tree.log')
    os.remove('tree.log')
    checkpoint = alnFile + '.ckp.gz'
    if os.path.exists(checkpoint):
        os.remove(checkpoint)

#***********************************************************************
def makeMidPointRootTree(treeFile):
    '''
     - Displays a dendogram of the tree generated from cluster representatives
    '''

    # Read the tree in newick format
    tree = Phylo.read(treeFile,'newick')

    # Root tree at midpoint
    tree.root_at_midpoint()

    # Write midpoint root tree
    Phylo.write(tree,'clsReps.aln.midpoint.treefile','newick')
    _validate_text_artifact('clsReps.aln.midpoint.treefile', "midpoint rooting")

#***********************************************************************

def drawAsciiTree(treeFile):
    '''
      - draws ASCII tree on the console
    '''

    if not os.path.exists(treeFile):
        raise ValueError("tree file <%s> does not exist" % treeFile)
    if os.path.getsize(treeFile) == 0:
        raise ValueError("tree file <%s> is empty" % treeFile)
    try:
        tree = Phylo.read(treeFile,'newick')
        Phylo.draw_ascii(tree)
    except Exception as exc:
        raise ValueError("could not draw Newick tree <%s>: %s" % (treeFile, exc)) from exc

#***********************************************************************
def alnLongSequenceClustersGNUParallel(nClusters,thread,mIterL):
    '''
        Align long sequence cluster files using Python-managed subprocesses.

        The legacy function name is retained for API compatibility.
    '''
    tasks = _long_alignment_tasks(nClusters)
    if not tasks:
        return

    def worker(task, threads_per_job):
        method, source, destination, record_ids, cluster = task
        if method == "copy":
            _copy_fasta(source, destination, expected_ids=record_ids)
        elif method == "linsi":
            linsi(source, destination, threads_per_job, mIterL, "align.log", cluster)
        else:
            fftnsi(source, destination, threads_per_job, mIterL, "align.log", cluster)

    _run_threaded(tasks, worker, thread)

#***********************************************************************
def alnLongSequenceClustersJoblibParallel(nClusters,thread,mIterL):
  '''
    - align long sequence cluster files in parallel using joblib
  '''

  log = 'align.log'
  #lh = open('pipelign.log','a')

  msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
  msg += ' Started aligning long sequence clusters\n'
  print(msg)

  tasks = _long_alignment_tasks(nClusters)
  if not tasks:
    return
  jobs, threads_per_job = _parallel_worker_plan(len(tasks), thread)

  def run_task(task):
    method, source, destination, record_ids, cluster = task
    if method == "copy":
      _copy_fasta(source, destination, expected_ids=record_ids)
    elif method == "linsi":
      linsi(source, destination, threads_per_job, mIterL, log, cluster)
    else:
      fftnsi(source, destination, threads_per_job, mIterL, log, cluster)

  joblib.Parallel(n_jobs=jobs, prefer="threads")(
    joblib.delayed(run_task)(task) for task in tasks
  )

  msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
  msg += ' Finished aligning long sequence clusters. Created files: <long.x.aln>\n'
  print(msg)


#***********************************************************************

def _long_alignment_tasks(nClusters):
  tasks = []
  for cluster in range(nClusters):
    source = 'long.' + str(cluster) + '.fas'
    destination = 'long.' + str(cluster) + '.aln'
    records = _validate_fasta_artifact(source, "cluster preparation")
    count = len(records)
    record_ids = [record.id for record in records]
    if _reusable_fasta(
      destination,
      alignment=True,
      expected_ids=record_ids,
      expected_ungapped=_record_sequences(records),
    ):
      continue
    method = "copy" if count == 1 else ("linsi" if count <= 100 else "fftnsi")
    tasks.append((method, source, destination, record_ids, str(cluster)))
  return tasks


#***********************************************************************


def makeHMMdbParallel(nClusters,thread,alphabet,log):
    '''
        - builds HMM from long sequence clusters
        - builds the database from HMMs
    '''

    # Fork the worker processes to perform computation concurrently
    # Create parameter map
    hmm = lambda i : ('long.' + str(i) + '.aln', 'long.' + str(i) + '.hmm')

    to_run_tuples = list(map(hmm, range(nClusters)))
    if _RESUME_EXISTING:
        to_run_tuples = [
            task for task in to_run_tuples
            if not _reusable_file(task[1])
        ]

    if len(to_run_tuples):
        msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Started building HMMs from long cluster alignments\n'
        print(msg)

        num_parallel_jobs, num_threads_per_job = _parallel_worker_plan(len(to_run_tuples), thread)
        joblib.Parallel(n_jobs=num_parallel_jobs, prefer="threads")(
            joblib.delayed(buildHMM)(
                x[0], x[1], num_threads_per_job, alphabet, log, x[0].split('.')[1]
            )
            for x in to_run_tuples
        )

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Building HMM database\n'
    print(msg)

    buildHMMdb(nClusters)

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' HMM database created\n'
    print(msg)

#***********************************************************************
def buildHMM(aName,hName,thread,alphabet,log,p=None):
    '''
        - Builds HMM for an alignment
    '''

    if p:
        lName = log + '.' + str(p)
    else:
        lName = log + '.temp'

    _validate_fasta_artifact(aName, "HMMER input", alignment=True)
    if _reusable_file(hName):
        try:
            _validate_hmm_artifact(hName, "resume")
            return
        except RuntimeError:
            pass

    gName = aName + '.temp'
    # create a temporary file to hold the alignment to build HMM
    shutil.copyfile(aName, gName)

    #*****testing for failure of hmmbuild
    # hmmbuild fails when single sequence is present
    # also when alignment does not have any gaps
    # solution: make a copy of the single sequnce
    # add a '-' at the end of each sequences

    # read in the alignment file
    aseqs = list(SeqIO.parse(gName,'fasta'))

    # for a single sequence cluster, add itself
    if len(aseqs) == 1:
        aseqs.append(copy.deepcopy(aseqs[0]))
        aseqs[1].id = aseqs[0].id + '_1'
        aseqs[1].name = aseqs[1].id
        aseqs[1].description = aseqs[1].id

    # for cluster with multiple sequences
    if len(aseqs) > 1:
        sumGap = 0 # for counting gaps

        for seq in aseqs:
            sumGap = sumGap + seq.seq.count('-')

        # if the alignment has no gaps, add one
        if sumGap == 0:
            for seq in aseqs:
                seq.seq = seq.seq + '-'

            # Update 'temp.aln' file
            os.remove(gName)
            SeqIO.write(aseqs,gName,'fasta')


    # create the command for hmmbuild
    if alphabet == 'dna':
        cl = ['hmmbuild', '--dna', '--cpu', str(thread), hName, gName]

    elif alphabet == 'rna':
        cl = ['hmmbuild', '--rna', '--cpu', str(thread), hName, gName]

    elif alphabet == 'aa':
        cl = ['hmmbuild', '--amino', '--cpu', str(thread), hName, gName]

    # run hmmbuild command
    try:
        with open(lName, 'a') as lh:
            subprocess.check_call(cl,stdout=lh,stderr=lh)
        _validate_hmm_artifact(hName, "hmmbuild")
    finally:
        if os.path.exists(gName):
            os.remove(gName)

    if not p:
        concatLogFiles('pipelign.log',lName)
        os.remove(lName)

#***********************************************************************

def buildHMMdb(nClusters):
    '''
        - builds HMM database from cluster HMMs
    '''

    dbName = 'pipelign.hmm'
    if _RESUME_EXISTING:
        paths = [dbName] + [dbName + suffix for suffix in ('.h3f', '.h3i', '.h3m', '.h3p')]
        if all(_reusable_file(path) for path in paths):
            return

    # write all HMMs in one file
    with open(dbName,'w') as fh:
        for i in range(nClusters):
            hName = 'long.' + str(i) + '.hmm'
            with open(hName,'r') as hh:
                shutil.copyfileobj(hh,fh)

    # create HMM database
    cl = ['hmmpress', '-f', dbName]

    # run command
    with open('buildhmm.log', 'a') as lh:
        subprocess.check_call(cl, stdout=lh, stderr=lh)
    _validate_hmm_artifact(dbName, "HMM database build")
    for suffix in ('.h3f', '.h3i', '.h3m', '.h3p'):
        _validate_text_artifact(dbName + suffix, "hmmpress")

    concatLogFiles('pipelign.log','buildhmm.log')
    os.remove('buildhmm.log')

#***********************************************************************
def longSeqClusters(tempDirPath,cDir,tDirName,tFileName):
    '''
    Copies the temporary working directory into current working directory
    The directory will contain long sequence cluster alignments, HMMs, log files
    '''

    msg = '\nPipelign has created long sequence cluster alignments and HMMs.\n'

    # save the pipelign run files into current directory

    # no path provided for temporary directory
    if tempDirPath is None:
        try:
            wName = cDir + '/' + zName
            shutil.copytree(tDirName,wName)
        except OSError as e:
            print(e)
            cZip(cDir,tFileName,zName)

    msg += '\nAlignment files and HMMs can be found in <%s>\n' % zName
    msg += '\tLong sequence alignments have names <long.xx.aln>\n'
    msg += '\tHMM file written in <long.xx.hmm>\n'
    msg += '\tHMM database written in <pipelign.hmm>\n'
    print(msg)

    print('\nThank you for using Pipelign.\n')
    return

#***********************************************************************

def addNameToInternalNodes(tree):
    # annotate names to internal nodes
    edge = 0
    for node in tree.traverse():
        if not _node_is_leaf(node):
            node.name = "NODE_%d" %edge
            edge += 1
    return tree

#************************************************************************

def changeLeafNamesToClusterNames(tree,tp):
    # annotate name of the leaf nodes to get cluster names
    for leaf in _iter_tree_leaves(tree):
        lName = tp + '.' + leaf.name.split('_')[-1]
        leaf.name = lName
        #print(leaf.name)
    return tree

#************************************************************************

def mergeCherries_gnu_parallel(nCherries,cStr,iNodes,fPre,thread,mIterL,log):
    '''
     Merge pairs concurrently using Python-managed subprocesses.

     The legacy function name is retained for API compatibility.
    '''
    msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Merging pairs of clusters\n'
    print(msg)

    cherries = cStr.split()
    if len(cherries) != nCherries * 2 or len(iNodes) != nCherries:
        raise ValueError("merge task metadata is inconsistent")
    tasks = [
        (cherries[index * 2], cherries[index * 2 + 1], iNodes[index], index)
        for index in range(nCherries)
    ]

    def worker(task, threads_per_job):
        first, second, output, index = task
        mergePair(
            first + '.aln',
            second + '.aln',
            output + '.aln',
            threads_per_job,
            mIterL,
            log,
            "%s.%d" % (fPre, index),
        )

    _run_threaded(tasks, worker, thread)

#************************************************************************

def mergePair(aln1,aln2,alnFile,thread,mIterLong,log,p=None):
    '''
      - merge pair of alignments using MAFFT's --merge method
    '''

    if p:
        lName = log + '.' + str(p)
    else:
        lName = log + '.temp'

    seq1 = _validate_fasta_artifact(aln1, "MAFFT merge input", alignment=True)
    seq2 = _validate_fasta_artifact(aln2, "MAFFT merge input", alignment=True)

    lSeq1 = len(seq1)
    lSeq2 = len(seq2)

    suffix = '.' + str(p) if p else '.temp'
    input_name = 'merge-input' + suffix + '.fas'
    table_name = 'merge-table' + suffix + '.txt'
    expected_ids = [record.id for record in seq1 + seq2]
    expected_ungapped = _record_sequences(seq1 + seq2)
    if _reusable_fasta(
        alnFile,
        alignment=True,
        expected_ids=expected_ids,
        expected_ungapped=expected_ungapped,
    ):
        return

    if lSeq1 == 1 and lSeq2 == 1:
        SeqIO.write(seq1 + seq2, input_name, 'fasta')
        cl = [
            'mafft', '--localpair', '--thread', str(thread), '--maxiterate',
            str(mIterLong), '--preservecase', input_name,
        ]
    elif lSeq1 == 1:
        cl = [
            'mafft', '--thread', str(thread), '--maxiterate', str(mIterLong),
            '--preservecase', '--addfragments', aln1, aln2,
        ]
    elif lSeq2 == 1:
        cl = [
            'mafft', '--thread', str(thread), '--maxiterate', str(mIterLong),
            '--preservecase', '--addfragments', aln2, aln1,
        ]
    else:
        SeqIO.write(seq1 + seq2, input_name, 'fasta')
        first_indices = ' '.join(str(index) for index in range(1, lSeq1 + 1))
        second_indices = ' '.join(
            str(index) for index in range(lSeq1 + 1, lSeq1 + lSeq2 + 1)
        )
        with open(table_name, 'w') as table:
            table.write(first_indices + '\n' + second_indices)
        cl = [
            'mafft', '--thread', str(thread), '--maxiterate', str(mIterLong),
            '--preservecase', '--merge', table_name, input_name,
        ]

    with open(lName, 'a') as lh, open(alnFile, 'w') as ah:
        subprocess.check_call(cl,stdout=ah,stderr=lh)
    _validate_fasta_artifact(
        alnFile, "MAFFT merge", alignment=True, expected_ids=expected_ids
        , expected_ungapped=expected_ungapped
    )
    for temporary in (input_name, table_name):
        if os.path.exists(temporary):
            os.remove(temporary)

    if not p:
        concatLogFiles('pipelign.log',lName)
        os.remove(lName)

#******************************************

def merge_gnu_parallel(repTree,tp,thread,mIterL,log):
    # read the tree
    tree = _read_named_tree(repTree)

    # annotate internal node names
    tree = addNameToInternalNodes(tree)

    # display the tree
    #print(tree.get_ascii(show_internal=True))

    # change leaf node names as cluster names
    if tp == 'long':
        tree = changeLeafNamesToClusterNames(tree,'long')
    elif tp == 'all':
        tree = changeLeafNamesToClusterNames(tree,'cls')
    print(tree)

    # get cherries and update trees
    while(len(tree) > 1):
        cherries = [] # empty list for cherries
        iNodes = [] # empty list for internal nodes
        for leaf in _iter_tree_leaves(tree):
            parent = getattr(leaf, "up", None)
            if parent is None:
                continue
            children = list(parent)
            if len(children) == 2 and all(_node_is_leaf(child) for child in children): # two leaf nodes
                #print(leaf.name)
                #print(leaf.up)
                cherries.append(leaf.name)
                if parent.name not in iNodes:
                    iNodes.append(parent.name)

        # call merge()
        if not cherries:
            raise RuntimeError("no mergeable cherries found in <%s>" % repTree)
        if len(cherries) % 2 != 0:
            raise PipelignError('incomplete cherries found: exiting')
        # even number of leafs in cherries list
        # get number of cherries
        nCherries = int(len(cherries) / 2)
        #print(nCherries)
        cStr = ' '.join(cherries)
        #print(cStr)
        mergeCherries_gnu_parallel(nCherries,cStr,iNodes,'merge',thread,mIterL,log)
        #break

        # get the pointers to the cherries leafs and remove them
        for cherry in cherries:
            l = _find_tree_node(tree, cherry)
            #print(l)
            #print(tree.get_ascii(show_internal=True))
            l.detach()
            #print(tree.get_ascii(show_internal=True))

    fName = tree.name + '.aln'
    return fName


#************************************************************************
def addLongSeqs(sName,aName,oName,thread,mIterL,log):
    '''
        - add long sequences to the alignment using MAFFT's --addfull
    '''

    lName = log + '.temp'

    sequence_records = _validate_fasta_artifact(sName, "MAFFT full-sequence input")
    alignment_records = _validate_fasta_artifact(aName, "MAFFT alignment input", alignment=True)
    expected_ids = [record.id for record in alignment_records + sequence_records]
    expected_ungapped = _record_sequences(alignment_records + sequence_records)
    if _reusable_fasta(
        oName,
        alignment=True,
        expected_ids=expected_ids,
        expected_ungapped=expected_ungapped,
    ):
        return

    cl = ['mafft', '--preservecase', '--thread', str(thread), '--maxiterate', str(mIterL), '--addfull', sName, aName]

    with open(lName, 'a') as lh, open(oName, 'w') as ah:
        subprocess.check_call(cl,stdout=ah,stderr=lh)
    _validate_fasta_artifact(
        oName,
        "MAFFT full-sequence alignment",
        alignment=True,
        expected_ids=expected_ids,
    )

    concatLogFiles(log,lName)
    os.remove(lName)

#***********************************************************************

def mergeConsensus(numClusters,oName,tp,mArgs,log):
    '''
        This function merges clusters together by following:
            - first consensus sequence is created from each cluster alignments
            - these consensus sequences are aligned together using MAFFT
            - gaps are placed in the original alignments based on the gaps in the consensus alignment

        tp can be:
            - 'long' for long sequence alignment
            - 'all' for cluster alignments with fragments

        Final alignment is produced by sequentially adding each of the clusters in the output
    '''

    ## get the file names prefix
    if tp == 'all':
        tp = 'cls'

    ## create a string to hold the consensus sequences in a FASTA format
    st = ''

    ## create consensus from each of the cluster alignments and create a sequence file
    for i in range(numClusters):
        aName = tp + '.' + str(i) + '.aln'
        with open(aName) as ah:
            aln = AlignIO.read(ah,'fasta')
        con = consensus_from_alignment(aln, mArgs.alphabet)

        st += '>' + aName + '\n'
        st += str(con) + '\n'

    ## write the consensus sequence file
    with open('consensus.fas','w') as fh:
        fh.write(st)

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Aligning consensus sequences from each clusters\n'
    print(msg)

    ## align the consensus sequences using FFT-NS-i
    fftnsi('consensus.fas','consensus.aln',mArgs.thread,mArgs.mIterL,log)

    ## read in the consensus alignment
    ## get positions for gaps in each sequence
    ## add gaps to their respective positions in the cluster alignments
    expected_ids = []
    with open(oName,'w') as fh: # create the output alignment file
        # read consensus alignment
        with open('consensus.aln') as ch:
            aln = AlignIO.read(ch,'fasta')

        for i in range(numClusters):
            indices = [j for j, x in enumerate(aln[i].seq) if x == "-"]
            cName = aln[i].id

            # read in the cluster alignment file
            with open(cName) as cHandle:
                seqs = SeqIO.parse(cHandle,'fasta')

                for seq in seqs:
                    expected_ids.append(seq.id)
                    # create a string object to write in the file
                    st = '>' + seq.id + '\n'

                    start = 0

                    for j in range(len(indices)):
                        end = indices[j] - j
                        st += seq.seq[start:end] + '-'
                        start = end

                    st += seq.seq[start:] + '\n'

                    fh.write(str(st))
    _validate_fasta_artifact(
        oName,
        "consensus merge",
        alignment=True,
        expected_ids=expected_ids,
    )


#***********************************************************************

def longSeqAlignmentConsensus(numClusters,mArgs):
    '''
    Creates alignment from only long sequence clusters
      - adds sequences from <long.bad.fas> if -b flag is used

    '''

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Merging long sequence clusters using consensus\n'
    print(msg)

    if numClusters == 1: # only one cluster present
        copyFile('long.0.aln','long.noBadSeq.aln')
    else:
        mergeConsensus(numClusters,'long.noBadSeq.aln','long',mArgs,'consensus.log')

    # check if <long.bad.fas> needs to be added to the final alignment
    if mArgs.keepBadSeqs and checkPresenceOfFile('long.bad.fas'):
        addLongSeqs('long.bad.fas','long.noBadSeq.aln','final.all.aln',mArgs.thread,mArgs.mIterL,'addBadSeq.log')
    else:
        copyFile('long.noBadSeq.aln','final.all.aln')

    '''
    elif numClusters == 2: # two clusters
        #mergePair('long.0.aln','long.1.aln','long.noBadSeq.aln',mArgs.thread,mArgs.mIterM,'merge.log',cDir,tFileName,zName)


    elif numClusters >= 3: # more than 2 clusters
        #resName = merge_gnu_parallel('clsReps.aln.midpoint.treefile','long',mArgs.thread,mArgs.mIterL,'merge.log',cDir,tFileName,zName)
        #copyFile(resName, 'long.noBadSeq.aln')
    '''
    '''
    # check if <long.bad.fas> needs to be added to the final alignment
    if mArgs.keepBadSeqs and checkPresenceOfFile('long.bad.fas'):
        addLongSeqs('long.bad.fas','long.noBadSeq.aln','final.aln',mArgs.thread,mArgs.mIterL,'addBadSeq.log',cDir,tFileName,zName)
    else:
        copyFile('long.noBadSeq.aln','final.aln')

    if checkPresenceOfFile('final.aln'):
        copyFile('final.aln',mArgs.outFile)


    if checkPresenceOfFile(mArgs.outFile):
        msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Final alignment written in <%s>\n' % mArgs.outFile
        print(msg)
    else:
        msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Pipelign could not write the file <%s>\n' % mArgs.outFile
        msg += '\tSomething went wrong. Please check the zipped output directory.\n'
        print(msg)

    '''
    # print('\nThank you for using Pipelign.\n')
    # sys.exit(0)

#***********************************************************************
def longSeqAlignmentParallel(numClusters,mArgs):
    '''
    Creates alignment from only long sequence clusters
      - adds sequences from <long.bad.fas> if -b flag is used

    '''

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Merging long sequence clusters using parallel strategy\n'
    print(msg)

    if numClusters == 1: # only one cluster present
        copyFile('long.0.aln','long.noBadSeq.aln')

    elif numClusters == 2: # two clusters
        mergePair('long.0.aln','long.1.aln','long.noBadSeq.aln',mArgs.thread,mArgs.mIterM,'merge.log')

    elif numClusters >= 3: # more than 2 clusters
        resName = merge_gnu_parallel('clsReps.aln.midpoint.treefile','long',mArgs.thread,mArgs.mIterL,'merge.log')
        copyFile(resName, 'long.noBadSeq.aln')

    # check if <long.bad.fas> needs to be added to the final alignment
    if mArgs.keepBadSeqs and checkPresenceOfFile('long.bad.fas'):
        addLongSeqs('long.bad.fas','long.noBadSeq.aln','final.all.aln',mArgs.thread,mArgs.mIterL,'addBadSeq.log')
    else:
        copyFile('long.noBadSeq.aln','final.all.aln')


#***********************************************************************

def runBlast(log,mArgs):
    '''
        - creates BLAST database from cluster representatives
        - search fragments against BLAST database
        - assign clusters to fragments and writes into file <frags.ClusterList.txt>
    '''
    # first create blast database of cluster reps
    bName = log + '.temp'
    cStr = 'makeblastdb -in clsReps.fas -input_type fasta -title pipelign -dbtype '
    if mArgs.alphabet in ['dna','rna']:
        cStr += 'nucl '
    elif mArgs.alphabet == 'aa':
        cStr += 'prot '

    cStr += '-out pipelign.blastdb'

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Creating BLAST database from Cluster representatives\n'
    print(msg)

    with open(bName, 'a') as bl:
        cl = cStr.split()
        subprocess.check_call(cl,stdout=bl,stderr=bl)
    database_files = [
        name for name in os.listdir('.')
        if name.startswith('pipelign.blastdb.')
    ]
    if not database_files:
        raise RuntimeError("makeblastdb did not produce database index files")

    # search fragments on the database to assign clusters
    cStr = ''
    if mArgs.alphabet in ['dna','rna']:
        cStr += 'blastn '
    elif mArgs.alphabet == 'aa':
        cStr += 'blastp '

    cStr += '-query frag.fas -db pipelign.blastdb -max_target_seqs 1 -outfmt 6 -evalue 0.05'
    #cStr += '-evalue 20 | sort -u -k1,2'

    with open('frag.blast.txt', 'w') as bo, open(bName, 'a') as bl:
        cl = cStr.split()
        subprocess.check_call(cl,stdout=bo,stderr=bl)
    if not os.path.isfile('frag.blast.txt'):
        raise RuntimeError("BLAST did not produce <frag.blast.txt>")

    # get uniq BLAST hits for each fragments



    # parse blast out file to create fragment clusterList file
    with open('frag.blast.txt','r') as fh:
        lines = [line.strip() for line in fh]

    st = '' # holds string for frag.ClusterList.txt

    mids = [] # holds names of the matched fragments
    matched = set()

    for line_number, line in enumerate(lines, 1):
        if not line:
            continue
        words = _fields(line, 2, "BLAST", line_number)
        cluster_number = words[1].rsplit('_', 1)[-1]
        if not cluster_number.isdigit():
            raise ValueError(
                "could not parse cluster number from BLAST subject %r" % words[1]
            )
        if words[0] not in matched:
            st += words[0] + '\t' + cluster_number + '\n'
            mids.append(words[0])
            matched.add(words[0])

    fh = open('frag.ClusterList.txt','w')
    fh.write(st)
    fh.close()

    concatLogFiles('pipelign.log',bName)
    os.remove(bName)

    return mids

#************************************************************************
def searchHMMdb(log,thread,alpha,res):
    '''
        HMM database is searched with the fragments to assign them a cliuster for alignment
    '''

    lName = log + '.temp'
    # generate the command for HMM search
    if alpha == 'dna' or alpha == 'rna':
        cStr = 'nhmmscan --cpu %d --tblout %s --noali  -E 0.05 pipelign.hmm frag.noBlast.fas' % (thread,res)

    elif alpha == 'aa':
        cStr = 'hmmscan --cpu %d --tblout %s --noali -E 0.05 pipelign.hmm frag.noBlast.fas' % (thread,res)

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Searching HMM database to assign clusters to remaining fragments\n'
    print(msg)

    with open(lName, 'a') as lh:
        cl = cStr.split()
        subprocess.check_call(cl,stdout=lh,stderr=lh)
    if not os.path.isfile(res):
        raise RuntimeError("HMMER did not produce <%s>" % res)

    concatLogFiles('pipelign.log',lName)
    os.remove(lName)

#**********************************************************************
def createFragmentClusters(numClusters):
    '''
        Fragments are written into separate files based on their cluster assignments
            - Reads in <frag.fas>
            - Reads in <frag.ClusterList.txt>
            - Creates fragments clusters <frag.x.fas>
            - Creates <frag.noClusters.fas> if no match found
    '''

    # read in the fragment sequnece file
    fseqs = list(SeqIO.parse('frag.fas','fasta'))

    fragments_by_id = _records_by_unique_id(fseqs, 'frag.fas')
    fids = list(fragments_by_id)
    fragment_index = {fragment_id: index for index, fragment_id in enumerate(fids)}

    # create an empty list to hold flags for fragments that are assigned clusters
    fragFlags = []
    for i in range(len(fids)):
        fragFlags.append(-1)


    # read in the cluster assignment file
    with open('frag.ClusterList.txt','r') as fh:
        clsList = [line.strip() for line in fh]

    ids = [] # holds IDs of the assigned fragments
    assigned_ids = set()
    cnum= [] # holds the cluster numbers

    for line_number, line in enumerate(clsList, 1):
        words = _fields(line, 2, 'frag.ClusterList.txt', line_number)
        if words[0] in assigned_ids:
            raise ValueError("duplicate fragment assignment for %r" % words[0])
        try:
            cluster_number = int(words[1])
        except ValueError as exc:
            raise ValueError("invalid fragment cluster number %r" % words[1]) from exc
        if cluster_number < 0 or cluster_number >= numClusters:
            raise ValueError("fragment cluster number %d is out of range" % cluster_number)
        ids.append(words[0])
        assigned_ids.add(words[0])
        cnum.append(cluster_number)
    #print(cnum)
    for i in range(numClusters):
        tfrags = [] # temporary list of all the fragments of that cluster

        # get the indices with cluster 'i'
        indices = [j for j, x in enumerate(cnum) if x == i]
        #print(indices)

        for ind in indices:
            fname = ids[ind] # get the sequence name/ID
            if fname in fragments_by_id: # if id matches to the fragment list
                sind = fragment_index[fname] # get the index in the fragment file
                tfrags.append(fragments_by_id[fname]) # add fragment to the cluster list
                fragFlags[sind] = 1 # flag

        if len(tfrags) > 0:
            fName = 'frag.' + str(i) + '.fas'
            SeqIO.write(tfrags,fName,'fasta')    # the the fragments to their cluster files

    # identify the orphan fragments and write them in a files

    orphans = []

    for i in range(len(fseqs)):
        if fragFlags[i] == -1:
            orphans.append(fseqs[i])

    if len(orphans) > 0:
        SeqIO.write(orphans,'frag.noClusters.fas','fasta')
        print("Unassigned fragments are written in <frag.noClusters.fas>\n")
        return 0

    else:
        return 1 # indicate that all fragments were assigned clusters

#************************************************************************
def getHMMclusters():
    '''
        Parse results from <hmm.out> to assign clusters to fragments based on HMM search
            - Reads in <frag.hmm.fas> file to read fragments searched against HMM_DB
            - Reads in <hmm.out> to parse output
            - Adds assigned fragments into the file <frag.ClusterList.txt>
            - Unassigned fragments are written into <frag.noClusters.fas>

    '''
    # open <frag.ClusterList.txt> for appending new entries
    # read in the hmm.out file
    with open('hmm.out','r') as hmm_handle:
        hmmOut = [line.strip() for line in hmm_handle]

    clusters_by_query = {}

    # make lists of queries and target HMMs
    for line_number, line in enumerate(hmmOut, 1):
        if not line.startswith('#'):
            if not line:
                continue
            words = _fields(line, 3, "HMMER", line_number)
            query = words[2]
            target = _cluster_number_from_name(words[0], "HMMER")
            if query not in clusters_by_query:
                clusters_by_query[query] = target

    # search each of the fragment ids in the qids list to find a match
    # get the first match as HMM search sorted the output based on best match
    with open('frag.ClusterList.txt','a') as cluster_handle:
        with open('frag.noBlast.fas','r') as handle:
            for seq in SeqIO.parse(handle,'fasta'):
                if seq.id in clusters_by_query: # if the id matched
                    cluster_handle.write("%s\t%s\n" %(seq.id, clusters_by_query[seq.id]))

    return len(clusters_by_query)

#**********************************************************************


def assignClusterstoFrags(numClusters,numFrag,mArgs):
    '''
    Assigns clusters to fragments based on BLAST and HMM search
    Copies the content of the temporary working directory into current working directory

    '''

    # if no fragments present in the dataset
    if not checkPresenceOfFile('frag.fas'):
        print('\nFragment file <frag.fas> does not exist. No fragments found for assigning clusters.\n')

        msg = '\nPipelign has created only long sequence cluster alignments and HMMs.\n'

    else: # <frag.fas> file exists, fragments need to be assigned clusters
        # running BLAST on the fragments
        mIds = runBlast('pipelign.log',mArgs) # running BLAST on the fragments
        fragWithClusters = len(mIds)
        #print(fragWithClusters)

        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' %d fragments were assigned clusters using BLAST search\n' % fragWithClusters
        print(msg)

        # not all the fragments were assigned clusters
        if fragWithClusters < int(numFrag):

            noMatchFrags = [] # hold fragments not matching in BLAST

            matched_ids = set(mIds)
            with open('frag.fas','r') as handle:
                for seq in SeqIO.parse(handle,'fasta'):
                    if seq.id not in matched_ids:
                        noMatchFrags.append(seq)

            if len(noMatchFrags) > 0:
                SeqIO.write(noMatchFrags,'frag.noBlast.fas','fasta')

            # search remaining fragments against HMM_DB
            searchHMMdb('pipelign.log', mArgs.thread, mArgs.alphabet, 'hmm.out')

            # get the clusters from HMM result
            numHMM = getHMMclusters()

            if numHMM > 0: # if fragments were assigned clusters based on HMM search
                msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
                msg += ' %d fragments were assigned clusters using HMM search\n' % numHMM
                print(msg)


        # create clusters of fragments, un-assigned ones are written in <frag.noCluster.fas>
        allAssigned = createFragmentClusters(numClusters)

#***********************************************************************

def addFragmentToClustersGNUParallel(nClusters,log,mArgs):
    '''
        Add fragments concurrently using Python-managed subprocesses.

        The legacy function name is retained for API compatibility.
    '''

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Adding fragments to long sequence clusters\n'
    print(msg)

    tasks = _fragment_alignment_tasks(nClusters)
    if not tasks:
        return

    def worker(task, threads_per_job):
        source, fragments, destination, cluster = task
        if fragments is None:
            records = _validate_fasta_artifact(source, "cluster alignment", alignment=True)
            _copy_fasta(
                source, destination, expected_ids=[record.id for record in records]
            )
        else:
            addFragments(
                fragments,
                source,
                destination,
                threads_per_job,
                mArgs.mIterL,
                log,
                cluster,
            )

    _run_threaded(tasks, worker, mArgs.thread)


#************************************************************************
def addFragments(fName,aName,oName,thread,mIterL,log,p=None):
    '''
        - add fragments to the alignment using MAFFT's --addfragments
    '''

    if p:
        lName = log + '.' + str(p)
    else:
        lName = log + '.temp'

    fragment_records = _validate_fasta_artifact(fName, "MAFFT fragment input")
    alignment_records = _validate_fasta_artifact(aName, "MAFFT alignment input", alignment=True)
    expected_ids = [record.id for record in alignment_records + fragment_records]
    expected_ungapped = _record_sequences(alignment_records + fragment_records)
    if _reusable_fasta(
        oName,
        alignment=True,
        expected_ids=expected_ids,
        expected_ungapped=expected_ungapped,
    ):
        return

    cl = ['mafft', '--preservecase', '--thread', str(thread), '--maxiterate', str(mIterL), '--addfragments', fName, aName]

    with open(lName, 'a') as lh, open(oName, 'w') as ah:
        subprocess.check_call(cl,stdout=ah,stderr=lh)
    _validate_fasta_artifact(
        oName,
        "MAFFT fragment alignment",
        alignment=True,
        expected_ids=expected_ids,
        expected_ungapped=expected_ungapped,
    )

    if not p:
        concatLogFiles('pipelign.log',lName)
        os.remove(lName)

#***********************************************************************

def addFragmentToClustersJoblibParallel(nClusters,log,thread,mIterL):
    '''
        This will add fragments to their cluster alignments using joblib
            - MAFFT's -addfragment will be used
            - will create cluster alignment files: <cls.x.aln>
    '''

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Adding fragments to long sequence clusters\n'
    print(msg)

    tasks = _fragment_alignment_tasks(nClusters)
    if not tasks:
        return
    jobs, threads_per_job = _parallel_worker_plan(len(tasks), thread)

    def run_task(task):
        source, fragments, destination, cluster = task
        if fragments is None:
            records = _validate_fasta_artifact(source, "cluster alignment", alignment=True)
            _copy_fasta(
                source, destination, expected_ids=[record.id for record in records]
            )
        else:
            addFragments(
                fragments, source, destination, threads_per_job, mIterL, log, cluster
            )

    joblib.Parallel(n_jobs=jobs, prefer="threads")(
        joblib.delayed(run_task)(task) for task in tasks
    )

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Finished adding fragments to long sequence clusters. Created files: <cls.x.aln>\n'
    print(msg)

#************************************************************************

def _fragment_alignment_tasks(nClusters):
    tasks = []
    for cluster in range(nClusters):
        source = 'long.' + str(cluster) + '.aln'
        fragment_path = 'frag.' + str(cluster) + '.fas'
        fragments = (
            fragment_path
            if os.path.isfile(fragment_path) and os.path.getsize(fragment_path) > 0
            else None
        )
        destination = 'cls.' + str(cluster) + '.aln'
        source_records = _validate_fasta_artifact(
            source, "cluster alignment", alignment=True
        )
        expected_ids = [record.id for record in source_records]
        expected_ungapped = _record_sequences(source_records)
        if fragments is not None:
            fragment_records = _validate_fasta_artifact(
                fragments, "fragment cluster"
            )
            expected_ids.extend(record.id for record in fragment_records)
            expected_ungapped.update(_record_sequences(fragment_records))
        if _reusable_fasta(
            destination,
            alignment=True,
            expected_ids=expected_ids,
            expected_ungapped=expected_ungapped,
        ):
            continue
        tasks.append((source, fragments, destination, str(cluster)))
    return tasks


#************************************************************************
def excludeClustersFromalignment(oAln,xAln):
    '''
    Excludes clusters from alignment
    '''
    # draw midpoint root tree
    print('\nMidpoint rooted tree for cluster representatives is drawn below. Tips are labelled as <H_N_C>.')
    print('\tH: header of the cluster representative sequence')
    print('\tN: number of long sequences in that cluster')
    print('\tC: unique identifier of the cluster')
    drawAsciiTree('clsReps.aln.midpoint.treefile')
    print('\n')


    # get user input for clusters
    inStr = input('Enter clusters (space separated) to exclude from final alignment:')
    inStr = inStr.strip().split()
    #print(inStr)

    if len(inStr) < 1:
        print('No cluster selected for exclusion')
        copyFile(oAln,xAln)
        return

    exList = list() # contains list of sequences for exclusion

    # read in the long cluster list file to list sequences for exclusion
    with open('long.ClusterList.txt','r') as fh:
        longSeqs = [line.strip() for line in fh]
    for lSeq in longSeqs:
        words = lSeq.split()
        if words[1] in inStr:
            exList.append(words[0])


    if checkPresenceOfFile('frag.ClusterList.txt'):
        # read in the fragment cluster list file to list fragments for exclusion
        with open('frag.ClusterList.txt','r') as fh:
            fragSeqs = [line.strip() for line in fh]
        for fSeq in fragSeqs:
            words = fSeq.split()
            if words[1] in inStr:
                exList.append(words[0])

    msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Removing sequences from final alignment\n'
    print(msg)

    if len(exList) == 0:
        print('No matching clusters selected for exclusion')
        copyFile(oAln, xAln)
        return

    if not checkPresenceOfFile(oAln):
        raise PipelignError('Original alignment file not found')
    aln = list(SeqIO.parse(oAln,'fasta'))

    excluded = set(exList)
    aSeq = [seq for seq in aln if seq.id not in excluded]
    if not aSeq:
        raise PipelignError('Cluster exclusion would remove every sequence')

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Writing new alignment file after removing sequences\n'
    print(msg)
    SeqIO.write(aSeq,xAln,'fasta')

#************************************************************************
# this is to check duplicate values for the same argument
class UniqueStore(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        if getattr(namespace, self.dest, self.default) is not None:
            print(option_string+values)
            parser.error(option_string + " appears several times.")
        setattr(namespace, self.dest, values)


#************************************************************************

def validateLengthThreshold(parser,lenThr):
    '''
    raises exception if length threshold is given multiple times
    if no value is given default value of 0.7 is returned

    '''

    if lenThr == None:
        print("== No value given for '-t', using the default value: 0.7")
        return 0.7
    elif len(lenThr) == 1:
        x = lenThr[0]
        if x < 0.5 or x > 1.0:
                parser.error('value given for -t is not in the range [0.5, 1.0]')
        return x
    else:
        parser.error('-t appears several times.')

#************************************************************************

def validateAlphabet(parser,alphabet):
    '''
    raises exception if alphabet is given multiple times
    if no value is given DNA is used by default

    '''
    if alphabet == None:
        print("== No value given for '-a', using the default sequence type: 'dna'")
        return 'dna'
    elif len(alphabet) == 1:
        return alphabet[0]
    else:
        parser.error('-a appears several times.')

#************************************************************************

def validateRunType(parser,run):
    '''
    raises exception if run type is given multiple times
    if no value is given, the built-in subprocess executor is used by default

    '''
    if run == None:
        print("== No value given for '-r', using built-in concurrent mode (G)")
        return 'G'
    elif len(run) == 1:
        #x = run[0].upper()
        #if x not in ['J','G']:
            #raise argparse.ArgumentTypeError('value of -r must be either J or G')
        #return x
        return run[0]

    else:
        parser.error('-r appears several times.')

#************************************************************************

def validateMergeType(parser,merge):
    '''
    raises exception if merge type is given multiple times
    if no value is given, merge is done parallelly by default

    '''
    if merge == None:
        print("== No value given for '-e', using default value: P")
        return 'P'
    elif len(merge) == 1:
        #x = run[0].upper()
        #if x not in ['J','G']:
            #raise argparse.ArgumentTypeError('value of -r must be either J or G')
        #return x
        return merge[0]

    else:
        parser.error('-r appears several times.')

#************************************************************************

def validatePercentageSimilarity(parser,simPer,alphabet):
    '''
    raises exception if percentage similarity threshold is given multiple times
    if no value is given default value of 0.8 is returned

    '''

    if simPer == None:
        print("== No value given for '-p', using the default value: 0.8")
        return 0.8
    elif len(simPer) == 1:
        x = float(simPer[0])
        if alphabet == 'aa':
            if x < 0.4 or x > 1.0:
                parser.error('value given for -p is not in the range [0.4, 1.0] for protein sequences')
            return x
        else:
            if x < 0.8 or x > 1.0:
                parser.error('value given for -p is not in the range [0.8, 1.0] for nucelotide sequences')
            return x

    else:
        parser.error('-p appears several times.')

#************************************************************************

def validateAmbigPercentage(parser,ambigPer):
    '''
    raises exception if percentage ambiguity threshold is given multiple times
    if no value is given default value of 0.1 is returned

    '''

    if ambigPer == None:
        print("== No value given for '-w', using the default value: 0.1")
        return 0.1
    elif len(ambigPer) == 1:
        x = float(ambigPer[0])
        if x < 0.0 or x > 1.0:
                parser.error('value given for -w is not in the range [0.0, 1.0]')
        return x
    else:
        parser.error('-p appears several times.')
#************************************************************************

def validateThread(parser,threads):
    '''
    checks number of available threads and sets value
    '''
    if threads == None:
        print("== No value given for '-q', running on a single CPU")
        return 1
    elif len(threads) == 1:
        if threads[0] < 1:
            parser.error('value given for -q must be a positive integer')
        numThreads = os.cpu_count()
        #print(numThreads,threads)
        if numThreads == None:
            return threads[0]
        elif numThreads >= threads[0]:
            #print(threads[0])
            return threads[0]
        elif numThreads < threads[0]:
            return numThreads
    else:
        parser.error('-q appears several times.')


#************************************************************************
def validateIterateMerge(parser,mergeIter):
    '''
    sets number of iterations for merging cluster alignments
    '''
    if mergeIter == None:
        print("== No value given for '-m', running only 1 iterations during merge")
        return 1
    elif len(mergeIter) == 1:
        if mergeIter[0] < 1:
            parser.error('value given for -m must be a positive integer')
        return mergeIter[0]
    else:
        parser.error('-m appears several times')

#************************************************************************
def validateIterateLong(parser,longIter):
    '''
    sets number of iterations for merging cluster alignments
    '''
    if longIter == None:
        print("== No value given for '-s', running only 1 iterations during aligning long sequence clusters")
        return 1
    elif len(longIter) == 1:
        if longIter[0] < 1:
            parser.error('value given for -s must be a positive integer')
        return longIter[0]
    else:
        parser.error('-s appears several times')

#************************************************************************
def getArguments(argv=None):
    '''
        Parses all the command line arguments from the user
    '''

    parser = argparse.ArgumentParser(description="Pipelign: creates multiple sequence alignment from FASTA formatted sequence file", formatter_class=argparse.RawTextHelpFormatter)
    #formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--inFile', required=True, help="Input sequence file in FASTA format",action=UniqueStore)
    parser.add_argument('-o', '--outFile', required=True, help="FASTA formatted output alignment file",action=UniqueStore)
    #parser.add_argument('-t', '--lenThr', type=lengthThreshold, help="Length threshold for full sequences (default: 0.9)",default=0.9,action='append')
    parser.add_argument('-t', '--lenThr', type=float, help="Length threshold for full sequences (default: 0.7)",action='append')
    #parser.add_argument('-c', '--code', type=int, help="Genetic code for translation (default: 1)",default=1, choices=[1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25])
    parser.add_argument('-a', '--alphabet', help='Input sequences can be dna/rna/aa (default: dna)', choices=['dna','aa','rna'], action='append')
    parser.add_argument('-f', '--keepOrphans', help='Add fragments without clusters', action="store_true")
    parser.add_argument('-b', '--keepBadSeqs', help='Add long sequences with too many ambiguous residues', action="store_true")
    parser.add_argument('-z', '--mZip', help='Create zipped intermediate output files', action="store_true")
    parser.add_argument('--resume', help='Resume a compatible run from its completed checkpoints', action="store_true")
    parser.add_argument('-p', '--simPer', type=float,help="Percent sequence similarity for clustering (default: 0.8)", action='append')
    parser.add_argument('-r', '--run', help="Run either (J)oblib or built-in concurrent mode (G, default)", choices=['J','G'],action='append')
    parser.add_argument('-e', '--merge', help="Merge using (P)arallel/(C)onsensus strategy  (default: P)", choices=['P','C'],action='append')
    parser.add_argument('-q', '--thread', type=int, help="Number of CPU/threads to use (default: 1)", action='append')
    parser.add_argument('-s', '--mIterateLong', type=int, help="Number of iterations to refine long alignments (default: 1)", action='append')
    parser.add_argument('-m', '--mIterateMerge', type=int, help="Number of iterations to refine merged alignment (default: 1)", action='append')
    parser.add_argument('-d', '--outDir', required=True, help="Name for output directory to hold intermediate files",action=UniqueStore)
    parser.add_argument('-c', '--clearExistingDirectory', help='Remove files from existing output directory', action="store_true")
    parser.add_argument('-w', '--ambigPer', type=float, help="Proportion of ambiguous characters allowed in the long sequences (default: 0.1)", action='append')
    parser.add_argument('-n', '--stage', type=int,default=5, choices=[1,2,3,4,5],
        help=textwrap.dedent('''\
        1  Make cluster alignments and HMM of long sequences
        2  Align long sequences only
        3  Assign fragments to clusters
        4  Make cluster alignments with fragments
        5  Align all sequences
  	    '''))

    parser.add_argument('-x', '--excludeClusters', help='Exclude clusters from final alignment', action="store_true")

    args = parser.parse_args(argv)

    # check arguments for validation

    args.lenThr = validateLengthThreshold(parser,args.lenThr)
    args.alphabet = validateAlphabet(parser,args.alphabet)
    args.simPer = validatePercentageSimilarity(parser,args.simPer,args.alphabet)
    args.run = validateRunType(parser,args.run)
    args.merge = validateMergeType(parser,args.merge)
    args.thread = validateThread(parser,args.thread)
    args.mIterateLong = validateIterateLong(parser,args.mIterateLong)
    args.mIterateMerge = validateIterateMerge(parser,args.mIterateMerge)
    args.ambigPer = validateAmbigPercentage(parser,args.ambigPer)

    #print(args.thread)
    #sys.exit(args.mIterateMerge)


    return args

#****************************************************************************
