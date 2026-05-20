from io import StringIO
import argparse

import pytest
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

import pipelign
from pipelign import doctor
from pipelign.cli import main as pipelign_main
from pipelign.gb2fas import convert_genbank_to_fasta
from pipelign.skill import SKILL_BYTES, skill_path


def _alignment(*seqs):
    fasta = "".join(">s%d\n%s\n" % (i, seq) for i, seq in enumerate(seqs))
    return AlignIO.read(StringIO(fasta), "fasta")


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
    assert ("parallel",) in stage_five
    assert ("blastn",) in stage_five


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
