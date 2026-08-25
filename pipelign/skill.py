"""Install a Copilot skill that teaches the pipelign golden path."""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass
from pathlib import Path


SKILL_NAME = "pipelign"
SKILL_FILENAME = "SKILL.md"
_COPILOT_HOME_ENV = "PIPELIGN_COPILOT_HOME"

SKILL_CONTENT = """---
name: pipelign
description: Use for pipelign workflows: automated multiple sequence alignment, FASTA/GenBank preprocessing, bioinformatics dependency checks, and troubleshooting MAFFT/CD-HIT/IQ-TREE/HMMER/BLAST installs.
argument-hint: "[alignment task or pipelign troubleshooting request]"
user-invocable: true
compatibility: Requires the pipelign CLI on PATH or a pipelign repository checkout.
---

# pipelign

Use this skill when the task involves `pipelign`, viral sequence alignment,
multiple sequence alignment pipelines, FASTA alignment, converting GenBank to
FASTA with `gb2fas`, or troubleshooting pipelign's external bioinformatics
tools.

Manual invocation is available as `/pipelign`.

## Golden path

1. Confirm which `pipelign` you are using.
   - Installed tool: `pipelign --help`
   - Repository checkout: `uvx --no-cache --from . pipelign --help`
2. Confirm the user's intended alphabet, run mode, and stage.
   - Alphabet: `dna`, `rna`, or `aa`
   - Run mode: `G` for built-in concurrency, `J` for Joblib
   - Stage: `1` through `5`
3. Check external executables with the same settings you will use for the run:
   ```bash
   pipelign-doctor -a dna -r G -n 5
   ```
4. If Python packaging is the problem, prefer a `uv` tool install for agents:
   ```bash
   uv tool install .
   ```
   For Git installs:
   ```bash
   uv tool install "pipelign @ git+https://github.com/asmmhossain/pipelign.git"
   ```
5. Install external bioinformatics binaries separately from conda-forge/bioconda:
   ```bash
   micromamba create -n pipelign -c conda-forge -c bioconda \
     python=3.12 mafft hmmer cd-hit iqtree blast
   micromamba activate pipelign
   ```
6. Start with a small input and explicit output paths:
   ```bash
   pipelign -i input.fasta -o aligned.fasta -d pipelign-work -a dna -r G -n 5
   ```
7. Report the final alignment path, intermediate directory, and any missing tools.

## Agent safety rules

- Do not assume the Python package install includes MAFFT, HMMER, CD-HIT,
  IQ-TREE, or BLAST. They are external executables.
- Run `pipelign-doctor` before invoking the full pipeline unless the user asks
  to skip checks.
- Use the user's intended `-a`, `-r`, and `-n` values in doctor checks.
- Prefer `-r G` for bounded built-in concurrency or `-r J` for Joblib.
- Keep experiments in an explicit working directory via `-d`.
- Do not overwrite an existing intermediate directory unless the user asks for
  `--clearExistingDirectory`.

## Useful commands

| Task | Command |
| --- | --- |
| Check all default DNA stage-5 dependencies | `pipelign-doctor -a dna -r G -n 5` |
| Check Joblib-only stage 1 dependencies | `pipelign-doctor -a dna -r J -n 1` |
| Convert GenBank to FASTA | `gb2fas input.gb output.fasta` |
| Run a DNA alignment | `pipelign -i input.fasta -o aligned.fasta -d work -a dna -r G -n 5` |
| Run with Joblib | `pipelign -i input.fasta -o aligned.fasta -d work -a dna -r J -n 5` |

## External tool packages

| Executable(s) | Conda package | Used for |
| --- | --- | --- |
| `mafft` | `mafft` | Sequence alignment and profile merging |
| `hmmbuild`, `hmmpress`, `hmmscan`, `nhmmscan` | `hmmer` | HMM building and sequence/profile search |
| `cd-hit`, `cd-hit-est` | `cd-hit` | Sequence clustering |
| `iqtree` or `iqtree2` | `iqtree` | Phylogenetic tree inference |
| `makeblastdb`, `blastn`, `blastp` | `blast` | Fragment-to-cluster assignment |

## Reporting back

When finished, report:

- The command path used and whether it was an installed tool or checkout.
- The doctor result for the user's intended settings.
- The output alignment path and intermediate work directory.
- Any missing external executables and the suggested conda package.
"""

SKILL_BYTES = SKILL_CONTENT.encode("utf-8")


@dataclass(frozen=True)
class SkillInstallResult:
    """Result of a skill install attempt."""

    status: str
    path: Path


class SkillInstallConflict(RuntimeError):
    """Raised when an existing user skill differs and --force was not used."""

    def __init__(self, path: Path) -> None:
        super().__init__(
            f"{path} already exists and differs from the pipelign skill. "
            "Use --force to overwrite it."
        )
        self.path = path


def default_copilot_home() -> Path:
    """Return the Copilot home directory used for personal skills."""

    configured = os.environ.get(_COPILOT_HOME_ENV)
    if configured:
        return Path(configured).expanduser()
    return Path.home() / ".copilot"


def skill_path(copilot_home: Path | None = None) -> Path:
    """Return the target SKILL.md path."""

    home = copilot_home.expanduser() if copilot_home is not None else default_copilot_home()
    return home / "skills" / SKILL_NAME / SKILL_FILENAME


def install_skill(
    *,
    copilot_home: Path | None = None,
    force: bool = False,
    dry_run: bool = False,
) -> SkillInstallResult:
    """Install the pipelign Copilot skill into a personal Copilot skills dir."""

    path = skill_path(copilot_home)
    exists = path.exists()
    if exists:
        existing = path.read_bytes()
        if existing == SKILL_BYTES:
            return SkillInstallResult("unchanged", path)
        if not force:
            if dry_run:
                return SkillInstallResult("would-conflict", path)
            raise SkillInstallConflict(path)
        if dry_run:
            return SkillInstallResult("would-overwrite", path)
        status = "updated"
    else:
        if dry_run:
            return SkillInstallResult("would-install", path)
        status = "installed"

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(SKILL_BYTES)
    return SkillInstallResult(status, path)


def skill_main(argv: list[str] | None = None) -> int:
    """Handle ``pipelign skill`` subcommands."""

    parser = argparse.ArgumentParser(
        prog="pipelign skill",
        description="Install and manage the pipelign Copilot skill.",
    )
    subparsers = parser.add_subparsers(dest="command")

    install_parser = subparsers.add_parser(
        "install",
        help="Install the pipelign personal Copilot skill.",
        description=(
            "Install SKILL.md into ~/.copilot/skills/pipelign/ so Copilot CLI "
            "can load pipelign golden-path guidance."
        ),
    )
    install_parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing pipelign skill if it differs.",
    )
    install_parser.add_argument(
        "--copilot-home",
        type=Path,
        default=None,
        help=(
            "Copilot home directory containing skills/ "
            f"(default: ${_COPILOT_HOME_ENV} or ~/.copilot)."
        ),
    )
    install_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be installed without writing files.",
    )

    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help(sys.stderr)
        return 2

    if args.command == "install":
        try:
            result = install_skill(
                copilot_home=args.copilot_home,
                force=args.force,
                dry_run=args.dry_run,
            )
        except SkillInstallConflict as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1
        if result.status == "unchanged":
            print(f"pipelign Copilot skill already up to date: {result.path}")
        elif result.status == "would-conflict":
            print(f"Would refuse to overwrite differing skill without --force: {result.path}")
        elif result.status == "would-overwrite":
            print(f"Would overwrite pipelign Copilot skill: {result.path}")
        elif result.status == "would-install":
            print(f"Would install pipelign Copilot skill: {result.path}")
        elif result.status == "updated":
            print(f"Updated pipelign Copilot skill: {result.path}")
        else:
            print(f"Installed pipelign Copilot skill: {result.path}")
        if result.status in {"installed", "updated", "unchanged"}:
            print("Reload Copilot CLI skills with: /skills reload")
            print("Inspect it with: /skills info pipelign")
        return 0

    parser.print_help(sys.stderr)
    return 2
