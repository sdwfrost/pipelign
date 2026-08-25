"""Check pipelign's external bioinformatics tool dependencies."""

import argparse
import sys

from pipelign import _first_available_executable, external_tool_details, required_external_tools


def main(argv=None):
    parser = argparse.ArgumentParser(description="Check external executables required by pipelign.")
    parser.add_argument("-a", "--alphabet", choices=["dna", "rna", "aa"], default="dna")
    parser.add_argument("-r", "--run", choices=["G", "J"], default="G", help="Parallel mode: built-in subprocess executor (G) or Joblib (J)")
    parser.add_argument("-n", "--stage", type=int, choices=[1, 2, 3, 4, 5], default=5)
    args = parser.parse_args(argv)

    print("Checking pipelign tools for alphabet=%s, run=%s, stage=%s" % (args.alphabet, args.run, args.stage))
    print("")

    missing = []
    for names in required_external_tools(args.alphabet, args.run, args.stage):
        details = external_tool_details(names)
        resolved = _first_available_executable(*names)
        if resolved:
            print("OK      %-24s %s" % (details["display"], resolved))
        else:
            print("MISSING %-24s package=%s; used for %s; %s" % (
                details["display"],
                details["package"],
                details["purpose"],
                details["stages"],
            ))
            missing.append(details)

    if missing:
        packages = sorted({details["package"] for details in missing if details["package"] != "unknown"})
        if packages:
            print("\nInstall missing tools from bioconda/conda-forge, for example:", file=sys.stderr)
            print("  micromamba install -c conda-forge -c bioconda " + " ".join(packages), file=sys.stderr)
        else:
            print("\nInstall missing tools from bioconda/conda-forge, then rerun pipelign.", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
