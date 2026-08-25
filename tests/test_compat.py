from io import StringIO
import argparse
import json
import os
import subprocess
from types import SimpleNamespace
import zipfile

import pytest
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO

import pipelign
from pipelign import doctor
from pipelign.cli import main as pipelign_main
from pipelign.gb2fas import convert_genbank_to_fasta
from pipelign.reporting import RunReport
from pipelign.skill import SKILL_BYTES, skill_path


def _alignment(*seqs):
    fasta = "".join(">s%d\n%s\n" % (i, seq) for i, seq in enumerate(seqs))
    return AlignIO.read(StringIO(fasta), "fasta")


def _write_executable(path, content):
    path.write_text("#!/usr/bin/env python3\n" + content)
    path.chmod(0o755)


def test_gen_code_limit_accepts_valid_code():
    assert pipelign.genCodeLimit("11") == 11


def test_gen_code_limit_rejects_invalid_code():
    with pytest.raises(argparse.ArgumentTypeError):
        pipelign.genCodeLimit("99")


def test_consensus_allows_lowercase_ties():
    assert str(pipelign.consensus_from_alignment(_alignment("A", "C"), "dna")) == "a"


def test_consensus_ignores_gaps_as_calls():
    assert str(pipelign.consensus_from_alignment(_alignment("-", "-", "A"), "dna")) == "A"


def test_consensus_rejects_empty_alignment():
    with pytest.raises(ValueError, match="empty alignment"):
        pipelign.consensus_from_alignment(MultipleSeqAlignment([]), "dna")


def test_consensus_all_gap_column_uses_ambiguous_call():
    assert str(pipelign.consensus_from_alignment(_alignment("-", "-"), "dna")) == "N"


def test_consensus_single_sequence_replaces_gap_with_call():
    assert str(pipelign.consensus_from_alignment(_alignment("AC-T"), "dna")) == "ACNT"


def test_consensus_protein_tie_uses_lowercase_call():
    assert str(pipelign.consensus_from_alignment(_alignment("A", "C"), "aa")) == "a"


def test_validate_fasta_type_accepts_rna(tmp_path):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">rna\nACGUACGU\n")

    assert pipelign.validateFastaType(fasta, "rna")
    assert not pipelign.validateFastaType(fasta, "dna")


def test_validate_fasta_type_rejects_dna_as_rna(tmp_path):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">dna\nACGTACGT\n")

    assert pipelign.validateFastaType(fasta, "dna")
    assert not pipelign.validateFastaType(fasta, "rna")


def test_dealign_accepts_single_sequence(tmp_path):
    input_path = tmp_path / "input.fasta"
    output_path = tmp_path / "output.fasta"
    input_path.write_text(">one\nAC-G.U\n")

    pipelign.deAlign(input_path, output_path)

    assert output_path.read_text() == ">one\nACGU\n"


def test_dealign_rejects_duplicate_identifiers(tmp_path):
    input_path = tmp_path / "input.fasta"
    output_path = tmp_path / "output.fasta"
    input_path.write_text(">duplicate\nACGU\n>duplicate\nACGU\n")

    with pytest.raises(ValueError, match="duplicate FASTA identifier"):
        pipelign.deAlign(input_path, output_path)


def test_separate_long_ambigs_exits_when_all_sequences_are_filtered(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "long.fas").write_text(">bad\nNNNN\n")

    with pytest.raises(pipelign.PipelignError, match="All long sequences"):
        pipelign.separateLongAmbigs("dna", 0.1)

    assert not (tmp_path / "long.good.fas").exists()


def test_ete4_tree_helpers_support_merge_operations(tmp_path):
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("((A:1,B:1)N1:1,(C:1,D:1)N2:1)Root;\n")

    tree = pipelign._read_named_tree(str(tree_file))
    leaf = pipelign._find_tree_node(tree, "A")

    assert leaf.name == "A"
    assert leaf.up.name == "N1"
    assert len(leaf.up) == 2

    leaf.detach()
    assert len(tree) == 3


def test_find_tree_node_rejects_duplicate_names(tmp_path):
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("((A:1,A:1)N1:1,B:1)Root;\n")
    tree = pipelign._read_named_tree(str(tree_file))

    with pytest.raises(ValueError, match="not unique"):
        pipelign._find_tree_node(tree, "A")


def test_draw_ascii_tree_rejects_empty_tree_file(tmp_path):
    tree_file = tmp_path / "empty.nwk"
    tree_file.write_text("")

    with pytest.raises(ValueError, match="empty"):
        pipelign.drawAsciiTree(str(tree_file))


def test_required_external_tools_are_stage_aware():
    stage_one = pipelign.required_external_tools("dna", "J", 1)
    stage_five = pipelign.required_external_tools("dna", "G", 5)

    assert ("parallel",) not in stage_one
    assert ("blastn",) not in stage_one
    assert ("parallel",) not in stage_five
    assert ("blastn",) in stage_five


@pytest.mark.parametrize(
    ("tasks", "threads", "expected"),
    [
        (100, 4, (4, 1)),
        (2, 64, (2, 32)),
        (5, 2, (2, 1)),
    ],
)
def test_parallel_worker_plan_respects_thread_limit(tasks, threads, expected):
    assert pipelign._parallel_worker_plan(tasks, threads) == expected


def test_builtin_alignment_worker_preserves_mafft_failure(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "long.0.fas").write_text(">a\nACGU\n>b\nACGU\n")
    error = subprocess.CalledProcessError(2, ["mafft"])
    monkeypatch.setattr(
        pipelign.subprocess,
        "check_call",
        lambda *args, **kwargs: (_ for _ in ()).throw(error),
    )

    with pytest.raises(subprocess.CalledProcessError) as exc_info:
        pipelign.alnLongSequenceClustersGNUParallel(1, 4, 1)

    assert exc_info.value is error
    assert not list(tmp_path.glob("*.sh"))


def test_build_hmm_duplicates_single_record_without_aliasing(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "single.aln").write_text(">original\nACGU\n")
    captured_ids = []

    def fake_hmmbuild(command, **kwargs):
        captured_ids.extend(
            record.id for record in SeqIO.parse(command[-1], "fasta")
        )
        (tmp_path / command[-2]).write_text("HMMER3/f\n")

    monkeypatch.setattr(pipelign.subprocess, "check_call", fake_hmmbuild)
    pipelign.buildHMM("single.aln", "single.hmm", 1, "rna", "hmmer.log", "0")

    assert captured_ids == ["original", "original_1"]


def test_mafft_zero_exit_without_output_is_rejected(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "input.fas").write_text(">one\nACGU\n>two\nACGU\n")
    monkeypatch.setattr(pipelign.subprocess, "check_call", lambda *args, **kwargs: None)

    with pytest.raises(RuntimeError, match="nonempty"):
        pipelign.linsi("input.fas", "output.aln", 1, 1, "mafft.log", "0")


def test_hmm_clusters_append_assignments(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "frag.ClusterList.txt").write_text("blast-hit\t0\n")
    (tmp_path / "frag.noBlast.fas").write_text(">hmm-hit\nACGU\n")
    (tmp_path / "hmm.out").write_text("long.2.aln - hmm-hit -\n")

    assert pipelign.getHMMclusters() == 1
    assert (tmp_path / "frag.ClusterList.txt").read_text() == (
        "blast-hit\t0\n"
        "hmm-hit\t2\n"
    )


def test_hmm_clusters_reject_malformed_output(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "frag.ClusterList.txt").write_text("")
    (tmp_path / "frag.noBlast.fas").write_text(">hmm-hit\nACGU\n")
    (tmp_path / "hmm.out").write_text("malformed\n")

    with pytest.raises(ValueError, match="malformed HMMER"):
        pipelign.getHMMclusters()


def test_blast_failure_propagates_original_process_error(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    error = subprocess.CalledProcessError(2, ["makeblastdb"])
    monkeypatch.setattr(pipelign.subprocess, "check_call", lambda *args, **kwargs: (_ for _ in ()).throw(error))

    with pytest.raises(subprocess.CalledProcessError) as exc_info:
        pipelign.runBlast("pipelign.log", SimpleNamespace(alphabet="dna"))

    assert exc_info.value is error


@pytest.mark.parametrize(
    "option",
    [
        ["-q", "0"],
        ["-q", "-1"],
        ["-s", "0"],
        ["-m", "-1"],
    ],
)
def test_cli_rejects_nonpositive_resource_values(option, tmp_path):
    argv = [
        "-i", str(tmp_path / "input.fasta"),
        "-o", str(tmp_path / "output.fasta"),
        "-d", str(tmp_path / "work"),
        *option,
    ]

    with pytest.raises(SystemExit) as exc_info:
        pipelign.getArguments(argv)

    assert exc_info.value.code == 2


def test_publish_alignment_is_atomic_and_validated(tmp_path):
    source = tmp_path / "source.aln"
    destination = tmp_path / "published.aln"
    source.write_text(">a\nAC-G\n>b\nACGG\n")

    pipelign.publishAlignment(source, destination)

    assert destination.read_text() == source.read_text()
    assert not (tmp_path / "published.aln.tmp").exists()
    assert not list(tmp_path.glob(".published.aln-*.tmp"))


def test_excluding_unknown_cluster_preserves_alignment(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "clsReps.aln.midpoint.treefile").write_text("(a:1,b:1);\n")
    (tmp_path / "long.ClusterList.txt").write_text("a\t0\nb\t1\n")
    (tmp_path / "original.aln").write_text(">a\nAC-G\n>b\nACGG\n")
    monkeypatch.setattr("builtins.input", lambda prompt: "99")

    pipelign.excludeClustersFromalignment("original.aln", "filtered.aln")

    assert (tmp_path / "filtered.aln").read_text() == (
        tmp_path / "original.aln"
    ).read_text()


def test_run_report_resumes_matching_input_and_archives(tmp_path, monkeypatch):
    input_path = tmp_path / "input.fas"
    output_path = tmp_path / "output.aln"
    work = tmp_path / "work"
    work.mkdir()
    input_path.write_text(">one\nACGU\n")
    arguments = {
        "alphabet": "rna",
        "lenThr": 0.7,
        "simPer": 0.8,
        "ambigPer": 0.1,
        "merge": "P",
        "resume": False,
    }
    monkeypatch.setattr("pipelign.reporting.tool_metadata", lambda tools: [])

    report = RunReport(work, input_path, output_path, arguments, [], archive=True)
    report.checkpoint("clustering", clusters=1)
    report.finish("success")

    data = json.loads((work / "pipelign-run.json").read_text())
    assert data["status"] == "success"
    assert data["checkpoints"]["clustering"]["clusters"] == 1
    assert (tmp_path / "work.zip").exists()
    with zipfile.ZipFile(tmp_path / "work.zip") as archive:
        archived_report = json.loads(archive.read("pipelign-run.json"))
    assert archived_report["archive"] == str(tmp_path / "work.zip")

    resumed_arguments = dict(arguments, resume=True)
    resumed = RunReport.resume(
        work, input_path, output_path, resumed_arguments, archive=False
    )
    assert resumed.has_checkpoint("clustering")
    resumed.finish("success")
    resumed_data = json.loads((work / "pipelign-run.json").read_text())
    assert resumed_data["archive"] is None


def test_run_report_rejects_changed_alignment_iterations(tmp_path, monkeypatch):
    input_path = tmp_path / "input.fas"
    output_path = tmp_path / "output.aln"
    work = tmp_path / "work"
    work.mkdir()
    input_path.write_text(">one\nACGU\n")
    arguments = {
        "alphabet": "rna",
        "lenThr": 0.7,
        "simPer": 0.8,
        "ambigPer": 0.1,
        "merge": "P",
        "mIterateLong": 1,
        "mIterateMerge": 1,
    }
    monkeypatch.setattr("pipelign.reporting.tool_metadata", lambda tools: [])
    RunReport(work, input_path, output_path, arguments, []).finish("failed")

    with pytest.raises(ValueError, match="mIterateLong"):
        RunReport.resume(
            work,
            input_path,
            output_path,
            dict(arguments, mIterateLong=2),
        )


def test_resume_rejects_same_id_alignment_with_changed_residues(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    alignment = tmp_path / "long.0.aln"
    alignment.write_text(">one\nAC-G\n")
    report_dir = tmp_path / "work"
    report_dir.mkdir()
    input_path = tmp_path / "input.fas"
    input_path.write_text(">one\nACG\n")
    monkeypatch.setattr("pipelign.reporting.tool_metadata", lambda tools: [])
    report = RunReport(report_dir, input_path, tmp_path / "out.aln", {}, [])
    (report_dir / "long.0.aln").write_text(alignment.read_text())
    report.checkpoint("alignment", artifacts=["long.0.aln"])
    checksums = report.artifact_checksums()
    report.finish("failed")

    (report_dir / "long.0.aln").write_text(">one\nATT\n")
    monkeypatch.chdir(report_dir)
    pipelign.configure_resume(True, checksums)
    try:
        assert not pipelign._reusable_fasta(
            "long.0.aln",
            alignment=True,
            expected_ids=["one"],
            expected_ungapped={"one": "ACG"},
        )
    finally:
        pipelign.configure_resume(False)


def test_archive_failure_marks_report_failed(tmp_path, monkeypatch):
    input_path = tmp_path / "input.fas"
    output_path = tmp_path / "output.aln"
    work = tmp_path / "work"
    work.mkdir()
    input_path.write_text(">one\nACGU\n")
    monkeypatch.setattr("pipelign.reporting.tool_metadata", lambda tools: [])
    report = RunReport(work, input_path, output_path, {}, [], archive=True)
    monkeypatch.setattr(
        report,
        "_publish_archive",
        lambda: (_ for _ in ()).throw(OSError("disk full")),
    )

    with pytest.raises(OSError, match="disk full"):
        report.finish("success")

    data = json.loads((work / "pipelign-run.json").read_text())
    assert data["status"] == "failed"
    assert data["archive"] is None


def test_fake_mafft_integration(tmp_path, monkeypatch):
    executable = tmp_path / "mafft"
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "import pathlib, sys\n"
        "path = pathlib.Path(sys.argv[-1])\n"
        "sys.stdout.write(path.read_text())\n"
    )
    executable.chmod(0o755)
    input_path = tmp_path / "input.fas"
    output_path = tmp_path / "output.aln"
    input_path.write_text(">a\nACGU\n>b\nACGU\n")
    monkeypatch.setenv("PATH", str(tmp_path) + os.pathsep + os.environ["PATH"])

    pipelign.linsi(
        input_path, output_path, 1, 1, str(tmp_path / "mafft.log"), "0"
    )

    assert [record.id for record in SeqIO.parse(output_path, "fasta")] == ["a", "b"]


def test_stage_one_with_fake_external_tools_writes_report_and_archive(
    tmp_path, monkeypatch
):
    tools = tmp_path / "tools"
    tools.mkdir()
    _write_executable(
        tools / "cd-hit-est",
        "import pathlib, shutil, sys\n"
        "if '--version' in sys.argv:\n"
        "    print('cd-hit-est fake 1.0'); raise SystemExit(0)\n"
        "source = pathlib.Path(sys.argv[sys.argv.index('-i') + 1])\n"
        "output = pathlib.Path(sys.argv[sys.argv.index('-o') + 1])\n"
        "shutil.copyfile(source, output)\n"
        "ids = [line[1:].split()[0] for line in source.read_text().splitlines() "
        "if line.startswith('>')]\n"
        "clusters = ''.join(f'>Cluster {i}\\n0 4nt, >{name}... *\\n' "
        "for i, name in enumerate(ids))\n"
        "output.with_suffix(output.suffix + '.clstr').write_text(clusters)\n",
    )
    _write_executable(
        tools / "hmmbuild",
        "import pathlib, sys\n"
        "if '--version' in sys.argv:\n"
        "    print('hmmbuild fake 1.0'); raise SystemExit(0)\n"
        "pathlib.Path(sys.argv[-2]).write_text('HMMER3/f\\n')\n",
    )
    _write_executable(
        tools / "hmmpress",
        "import pathlib, sys\n"
        "if '--version' in sys.argv:\n"
        "    print('hmmpress fake 1.0'); raise SystemExit(0)\n"
        "database = pathlib.Path(sys.argv[-1])\n"
        "for suffix in ('.h3f', '.h3i', '.h3m', '.h3p'):\n"
        "    pathlib.Path(str(database) + suffix).write_text('index\\n')\n",
    )
    _write_executable(
        tools / "mafft",
        "import pathlib, sys\n"
        "if '--version' in sys.argv:\n"
        "    print('mafft fake 1.0'); raise SystemExit(0)\n"
        "paths = [pathlib.Path(sys.argv[-1])]\n"
        "if '--addfragments' in sys.argv:\n"
        "    paths.insert(0, pathlib.Path(sys.argv[sys.argv.index('--addfragments') + 1]))\n"
        "records = []\n"
        "for path in paths:\n"
        "    current = None\n"
        "    for line in path.read_text().splitlines():\n"
        "        if line.startswith('>'):\n"
        "            records.append([line, '']); current = records[-1]\n"
        "        elif current is not None:\n"
        "            current[1] += line\n"
        "width = max(len(sequence) for _, sequence in records)\n"
        "for header, sequence in records:\n"
        "    print(header); print(sequence.ljust(width, '-'))\n",
    )
    _write_executable(
        tools / "iqtree",
        "import pathlib, sys\n"
        "if '--version' in sys.argv:\n"
        "    print('iqtree fake 1.0'); raise SystemExit(0)\n"
        "alignment = pathlib.Path(sys.argv[sys.argv.index('-s') + 1])\n"
        "ids = [line[1:].split()[0] for line in alignment.read_text().splitlines() "
        "if line.startswith('>')]\n"
        "tree = f'(({ids[0]}:1,{ids[1]}:1):1,{ids[2]}:1);\\n'\n"
        "pathlib.Path(str(alignment) + '.treefile').write_text(tree)\n"
        "pathlib.Path(str(alignment) + '.ckp.gz').write_text('checkpoint')\n",
    )
    input_path = tmp_path / "input.fas"
    output_path = tmp_path / "unused.aln"
    work = tmp_path / "work"
    input_path.write_text(">alpha\nACGT\n>beta\nACGA\n")
    monkeypatch.setenv("PATH", str(tools) + os.pathsep + os.environ["PATH"])

    with pytest.raises(SystemExit) as exc_info:
        pipelign_main(
            [
                "-i", str(input_path),
                "-o", str(output_path),
                "-d", str(work),
                "-a", "dna",
                "-r", "G",
                "-n", "1",
                "-z",
            ]
        )

    assert exc_info.value.code == 0
    report = json.loads((work / "pipelign-run.json").read_text())
    assert report["status"] == "success"
    assert report["checkpoints"]["clustering"]["clusters"] == 2
    assert report["checkpoints"]["cluster_alignments_and_hmms"]["clusters"] == 2
    assert (tmp_path / "work.zip").exists()

    _write_executable(
        tools / "hmmbuild",
        "import sys\n"
        "if '--version' in sys.argv:\n"
        "    print('hmmbuild fake 1.0'); raise SystemExit(0)\n"
        "raise SystemExit(9)\n",
    )
    with pytest.raises(SystemExit) as resume_exc:
        pipelign_main(
            [
                "-i", str(input_path),
                "-o", str(output_path),
                "-d", str(work),
                "-a", "dna",
                "-r", "G",
                "-q", "2",
                "-n", "1",
                "--resume",
            ]
        )
    assert resume_exc.value.code == 0

    _write_executable(
        tools / "hmmbuild",
        "import pathlib, sys\n"
        "if '--version' in sys.argv:\n"
        "    print('hmmbuild fake 1.0'); raise SystemExit(0)\n"
        "pathlib.Path(sys.argv[-2]).write_text('HMMER3/f\\n')\n",
    )
    _write_executable(
        tools / "makeblastdb",
        "import pathlib, sys\n"
        "if '--version' in sys.argv:\n"
        "    print('makeblastdb fake 1.0'); raise SystemExit(0)\n"
        "source = pathlib.Path(sys.argv[sys.argv.index('-in') + 1])\n"
        "prefix = pathlib.Path(sys.argv[sys.argv.index('-out') + 1])\n"
        "subject = next(line[1:].split()[0] for line in source.read_text().splitlines() "
        "if line.startswith('>'))\n"
        "pathlib.Path(str(prefix) + '.nhr').write_text('index\\n')\n"
        "pathlib.Path(str(prefix) + '.subject').write_text(subject)\n",
    )
    _write_executable(
        tools / "blastn",
        "import pathlib, sys\n"
        "if '--version' in sys.argv:\n"
        "    print('blastn fake 1.0'); raise SystemExit(0)\n"
        "query = pathlib.Path(sys.argv[sys.argv.index('-query') + 1])\n"
        "prefix = pathlib.Path(sys.argv[sys.argv.index('-db') + 1])\n"
        "subject = pathlib.Path(str(prefix) + '.subject').read_text().strip()\n"
        "for line in query.read_text().splitlines():\n"
        "    if line.startswith('>'):\n"
        "        print(f'{line[1:].split()[0]}\\t{subject}')\n",
    )
    _write_executable(
        tools / "nhmmscan",
        "import sys\n"
        "if '--version' in sys.argv:\n"
        "    print('nhmmscan fake 1.0'); raise SystemExit(0)\n"
        "raise SystemExit(0)\n",
    )
    full_input = tmp_path / "full-input.fas"
    full_output = tmp_path / "full-output.aln"
    full_work = tmp_path / "full-work"
    full_input.write_text(
        ">alpha\nACGT\n>beta\nACGA\n>gamma\nACGC\n>fragment\nA\n"
    )

    with pytest.raises(SystemExit) as full_exc:
        pipelign_main(
            [
                "-i", str(full_input),
                "-o", str(full_output),
                "-d", str(full_work),
                "-a", "dna",
                "-r", "G",
                "-q", "2",
                "-n", "5",
            ]
        )

    assert full_exc.value.code == 0
    final_records = list(SeqIO.parse(full_output, "fasta"))
    assert {record.id for record in final_records} == {
        "alpha",
        "beta",
        "gamma",
        "fragment",
    }
    full_report = json.loads((full_work / "pipelign-run.json").read_text())
    assert full_report["status"] == "success"
    assert full_report["summary"]["final_sequences"] == 4
    assert full_report["checkpoints"]["final_alignment"]["output"] == str(full_output)

    joblib_output = tmp_path / "joblib-output.aln"
    joblib_work = tmp_path / "joblib-work"
    with pytest.raises(SystemExit) as joblib_exc:
        pipelign_main(
            [
                "-i", str(full_input),
                "-o", str(joblib_output),
                "-d", str(joblib_work),
                "-a", "dna",
                "-r", "J",
                "-q", "2",
                "-n", "5",
            ]
        )
    assert joblib_exc.value.code == 0
    assert {
        record.id for record in SeqIO.parse(joblib_output, "fasta")
    } == {record.id for record in final_records}


def test_gb2fas_rejects_empty_input(tmp_path):
    input_path = tmp_path / "empty.gb"
    output_path = tmp_path / "out.fas"
    input_path.write_text("")

    with pytest.raises(ValueError):
        convert_genbank_to_fasta(input_path, output_path)


def test_gb2fas_missing_output_directory_is_os_error(tmp_path):
    input_path = tmp_path / "input.gb"
    output_path = tmp_path / "missing" / "out.fas"
    input_path.write_text(
        "LOCUS       TEST0001       4 bp    DNA     linear   UNA 01-JAN-2000\n"
        "DEFINITION  Test record.\n"
        "ACCESSION   TEST0001\n"
        "VERSION     TEST0001.1\n"
        "FEATURES             Location/Qualifiers\n"
        "ORIGIN\n"
        "        1 acgt\n"
        "//\n"
    )

    with pytest.raises(OSError):
        convert_genbank_to_fasta(input_path, output_path)


def test_doctor_reports_actionable_missing_tools(monkeypatch, capsys):
    monkeypatch.setattr(doctor, "_first_available_executable", lambda *names: None)

    assert doctor.main(["-a", "dna", "-r", "J", "-n", "1"]) == 1

    out = capsys.readouterr().out
    assert "alphabet=dna, run=J, stage=1" in out
    assert "package=mafft" in out
    assert "parallel" not in out
    assert "blastn" not in out


def test_pipelign_cli_help_exits_zero(capsys):
    with pytest.raises(SystemExit) as exc_info:
        pipelign_main(["--help"])
    assert exc_info.value.code == 0
    assert "Pipelign: creates multiple sequence alignment" in capsys.readouterr().out


def test_pipelign_skill_install_writes_personal_skill(tmp_path, monkeypatch, capsys):
    copilot_home = tmp_path / "copilot"
    monkeypatch.setenv("PIPELIGN_COPILOT_HOME", str(copilot_home))

    assert pipelign_main(["skill", "install"]) == 0

    path = skill_path(copilot_home)
    assert path.exists()
    assert path.read_bytes() == SKILL_BYTES
    out = capsys.readouterr().out
    assert "Installed pipelign Copilot skill" in out
    assert "/skills reload" in out


def test_pipelign_skill_install_refuses_changed_file(tmp_path, capsys):
    copilot_home = tmp_path / "copilot"
    path = skill_path(copilot_home)
    path.parent.mkdir(parents=True)
    path.write_bytes(b"custom user skill\n")

    assert pipelign_main(["skill", "install", "--copilot-home", str(copilot_home)]) == 1

    err = capsys.readouterr().err
    assert "already exists and differs" in err
    assert path.read_bytes() == b"custom user skill\n"


def test_skill_subcommand_help(capsys):
    with pytest.raises(SystemExit) as exc_info:
        pipelign_main(["skill", "install", "--help"])
    assert exc_info.value.code == 0
    assert "--copilot-home" in capsys.readouterr().out
