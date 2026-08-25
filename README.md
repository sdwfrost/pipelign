# pipelign

A pipeline for automated multiple sequence alignment, particularly of viral
sequences.

Pipelign accepts FASTA input, separates full-length sequences from fragments,
clusters related sequences, builds cluster alignments and HMMs, assigns
fragments back to clusters, and produces a final multiple sequence alignment.

## Supported Python versions

The Python package supports Python 3.10 through 3.14. Runtime bioinformatics
executables are installed separately from conda-forge/bioconda or your system
package manager.

## External dependencies

Pipelign shells out to these tools:

| Executable(s) | Conda package | Used for |
| --- | --- | --- |
| `mafft` | `mafft` | Sequence alignment and profile merging |
| `hmmbuild`, `hmmpress`, `hmmscan`, `nhmmscan` | `hmmer` | HMM building and sequence/profile search |
| `cd-hit`, `cd-hit-est` | `cd-hit` | Sequence clustering |
| `iqtree` or `iqtree2` | `iqtree` | Tree inference |
| `makeblastdb`, `blastn`, `blastp` | `blast` | Fragment-to-cluster assignment |

Create an environment with micromamba:

```bash
micromamba create -n pipelign -c conda-forge -c bioconda \
  python=3.12 mafft hmmer cd-hit iqtree blast
micromamba activate pipelign
```

Or with conda:

```bash
conda create -n pipelign -c conda-forge -c bioconda \
  python=3.12 mafft hmmer cd-hit iqtree blast
conda activate pipelign
```

Or with pixi:

```bash
pixi init pipelign-env
cd pipelign-env
pixi add -c conda-forge -c bioconda python=3.12 mafft hmmer cd-hit iqtree blast
pixi shell
```

## Install the Python package

For agent-friendly use, install the command as a uv tool from a checkout:

```bash
git clone https://github.com/asmmhossain/pipelign.git
cd pipelign
uv tool install .
```

Update or remove the tool with:

```bash
uv tool upgrade pipelign
uv tool uninstall pipelign
```

For editable development inside an activated conda/micromamba/pixi environment:

```bash
uv pip install -e ".[dev]"
```

Plain pip also works:

```bash
python -m pip install .
python -m pip install -e ".[dev]"
```

## Copilot agent skill

Install the optional personal Copilot skill so agents get the recommended
pipelign workflow when relevant:

```bash
pipelign skill install
```

If Copilot CLI is already running, reload skills with `/skills reload` and
check it with `/skills info pipelign`.

## Check an installation

Run `pipelign-doctor` with the same alphabet, run mode, and stage you plan to
use:

```bash
pipelign-doctor -a dna -r G -n 5
```

Use `-r J` to check Joblib mode:

```bash
pipelign-doctor -a dna -r J -n 5
```

## Usage

Basic DNA alignment:

```bash
pipelign -i input.fasta -o aligned.fasta -d pipelign-work -a dna -r G -n 5
```

Run with Joblib:

```bash
pipelign -i input.fasta -o aligned.fasta -d pipelign-work -a dna -r J -n 5
```

Convert GenBank records to FASTA:

```bash
gb2fas input.gb output.fasta
```

### Main command reference

```text
usage: pipelign [-h] -i INFILE -o OUTFILE [-t LENTHR] [-a {dna,aa,rna}] [-f]
                [-b] [-z] [-p SIMPER] [-r {J,G}] [-e {P,C}] [-q THREAD]
                [-s MITERATELONG] [-m MITERATEMERGE] -d OUTDIR [-c]
                [-w AMBIGPER] [-n {1,2,3,4,5}] [-x]

Pipelign: creates multiple sequence alignment from FASTA formatted sequence file

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --inFile INFILE
                        Input sequence file in FASTA format
  -o OUTFILE, --outFile OUTFILE
                        FASTA formatted output alignment file
  -t LENTHR, --lenThr LENTHR
                        Length threshold for full sequences (default: 0.7)
  -a {dna,aa,rna}, --alphabet {dna,aa,rna}
                        Input sequences can be dna/rna/aa (default: dna)
  -f, --keepOrphans     Add fragments without clusters
  -b, --keepBadSeqs     Add long sequences with too many ambiguous residues
  -z, --mZip            Create zipped intermediate output files
  -p SIMPER, --simPer SIMPER
                        Percent sequence similarity for clustering (default: 0.8)
  -r {J,G}, --run {J,G}
                        Run either (J)oblib or built-in concurrent mode (G, default)
  -e {P,C}, --merge {P,C}
                        Merge using (P)arallel/(C)onsensus strategy (default: P)
  -q THREAD, --thread THREAD
                        Number of CPU/threads to use (default: 1)
  -s MITERATELONG, --mIterateLong MITERATELONG
                        Number of iterations to refine long alignments (default: 1)
  -m MITERATEMERGE, --mIterateMerge MITERATEMERGE
                        Number of iterations to refine merged alignment (default: 1)
  -d OUTDIR, --outDir OUTDIR
                        Name for output directory to hold intermediate files
  -c, --clearExistingDirectory
                        Remove files from existing output directory
  -w AMBIGPER, --ambigPer AMBIGPER
                        Proportion of ambiguous characters allowed in long sequences
  -n {1,2,3,4,5}, --stage {1,2,3,4,5}
                        1  Make cluster alignments and HMM of long sequences
                        2  Align long sequences only
                        3  Assign fragments to clusters
                        4  Make cluster alignments with fragments
                        5  Align all sequences
  -x, --excludeClusters
                        Exclude clusters from final alignment
```

## Changes from the original package

- Modern Python packaging via `pyproject.toml`.
- Python 3.10-3.14 compatibility.
- `ete4` replaces `ete3`.
- Biopython consensus generation uses `Bio.motifs` instead of deprecated
  `dumb_consensus()`.
- `pipelign`, `gb2fas`, and `pipelign-doctor` are installed as console scripts.
- `pipelign-doctor` reports missing external executables and install packages.
- `pipelign skill install` installs optional Copilot CLI workflow guidance.
- Built-in mode uses bounded Python-managed subprocess concurrency and does not
  require GNU Parallel.
- Every run writes `pipelign-run.json`; use `--resume` to reuse validated
  artifacts from a compatible interrupted run.
- `-z` writes an archive next to the intermediate directory after success.

## Development

Install development dependencies:

```bash
uv pip install -e ".[dev]"
```

Run tests and syntax checks:

```bash
make lint
make test
```

Build and check a distribution:

```bash
python -m build
python -m twine check dist/*
```

## Citation

Pipelign: an alignment pipeline for viral sequences. A.S.Md.M. Hossain and
S.D.W. Frost, in preparation.

## License

This Python package is licensed under the MIT license.

## Warranty

This software is provided "as is" and without any express or implied warranties,
including, without limitation, the implied warranties of merchantability and
fitness for a particular purpose.
