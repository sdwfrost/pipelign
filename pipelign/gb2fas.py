"""Convert GenBank files to accession-header FASTA."""

import argparse
import os
from pathlib import Path
import warnings

from Bio import BiopythonParserWarning, SeqIO


def convert_genbank_to_fasta(input_path, output_path):
    input_path = Path(input_path)
    output_path = Path(output_path)
    tmp_path = output_path.with_name(output_path.name + ".tmp")
    count = 0

    try:
        with input_path.open() as infile, tmp_path.open("w") as outfile:
            with warnings.catch_warnings():
                warnings.simplefilter("error", BiopythonParserWarning)
                for record in SeqIO.parse(infile, "genbank"):
                    record.id = record.name
                    record.description = record.name
                    outfile.write(record.format("fasta"))
                    count += 1
    except OSError:
        if tmp_path.exists():
            tmp_path.unlink()
        raise
    except Exception as exc:
        if tmp_path.exists():
            tmp_path.unlink()
        raise ValueError("could not parse GenBank file %s: %s" % (input_path, exc)) from exc

    if count == 0:
        if tmp_path.exists():
            tmp_path.unlink()
        raise ValueError("no GenBank records found in %s" % input_path)

    os.replace(tmp_path, output_path)
    return count


def main(argv=None):
    parser = argparse.ArgumentParser(description="Convert GenBank to accession-header FASTA.")
    parser.add_argument("input", help="Input GenBank file")
    parser.add_argument("output", help="Output FASTA file")
    args = parser.parse_args(argv)

    try:
        count = convert_genbank_to_fasta(args.input, args.output)
    except (OSError, ValueError) as exc:
        parser.exit(1, "gb2fas: error: %s\n" % exc)

    print("Converted %d record(s) from %s to %s" % (count, args.input, args.output))


if __name__ == "__main__":
    main()
