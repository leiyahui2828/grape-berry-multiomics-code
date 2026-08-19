#!/usr/bin/env python3

import argparse
import subprocess
import sys
import re
from concurrent.futures import ThreadPoolExecutor


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate sequence length, GC content and NMFE from FASTA using RNAfold."
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Input FASTA file"
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV file"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of parallel threads [default: 1]"
    )

    parser.add_argument(
        "--rnafold",
        type=str,
        default="RNAfold",
        help="Path to RNAfold executable [default: RNAfold]"
    )

    return parser.parse_args()


def detect_rnafold_noPS_arg(rnafold_exe):
    """
    Automatically detect whether RNAfold supports -noPS or --noPS.
    Old RNAfold usually uses -noPS.
    New RNAfold may support --noPS.
    """
    test_seq = "AUGGCUACGUAGCUAGCUAGC\n"

    for noPS_arg in ["-noPS", "--noPS"]:
        try:
            result = subprocess.run(
                [rnafold_exe, noPS_arg],
                input=test_seq,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
        except FileNotFoundError:
            sys.exit(
                f"Error: RNAfold was not found: {rnafold_exe}\n"
                "Please provide a valid RNAfold path using --rnafold."
            )

        if result.returncode == 0 and "(" in result.stdout and ")" in result.stdout:
            return noPS_arg

    sys.stderr.write("Error: failed to automatically detect RNAfold noPS parameter.\n")
    sys.stderr.write("Please test RNAfold manually:\n")
    sys.stderr.write(f"  echo 'AUGGCUACGUAGCUAGCUAGC' | {rnafold_exe} -noPS\n")
    sys.stderr.write(f"  echo 'AUGGCUACGUAGCUAGCUAGC' | {rnafold_exe} --noPS\n")
    sys.exit("RNAfold exists but cannot run correctly with -noPS or --noPS.")


def check_rnafold(rnafold_exe, noPS_arg):
    test_seq = "AUGGCUACGUAGCUAGCUAGC\n"

    try:
        result = subprocess.run(
            [rnafold_exe, noPS_arg],
            input=test_seq,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
    except FileNotFoundError:
        sys.exit(
            f"Error: RNAfold was not found: {rnafold_exe}\n"
            "Please provide a valid RNAfold path using --rnafold."
        )

    if result.returncode != 0:
        sys.stderr.write("RNAfold test command failed.\n")
        sys.stderr.write("Command: " + " ".join([rnafold_exe, noPS_arg]) + "\n")
        sys.stderr.write("STDOUT:\n" + result.stdout + "\n")
        sys.stderr.write("STDERR:\n" + result.stderr + "\n")
        sys.exit("Error: RNAfold exists but cannot run correctly.")

    if "(" not in result.stdout or ")" not in result.stdout:
        sys.stderr.write("Unexpected RNAfold test output:\n")
        sys.stderr.write(result.stdout + "\n")
        sys.exit("Error: RNAfold output format is unexpected.")


def read_fasta(fasta_file):
    seq_id = None
    seq_chunks = []

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(seq_chunks)

                seq_id = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if seq_id is not None:
            yield seq_id, "".join(seq_chunks)


def calc_gc_content(seq):
    seq = seq.upper()

    a = seq.count("A")
    t = seq.count("T")
    u = seq.count("U")
    g = seq.count("G")
    c = seq.count("C")

    valid_len = a + t + u + g + c

    if valid_len == 0:
        return None

    return (g + c) / valid_len


def parse_mfe_from_rnafold_output(rnafold_stdout):
    """
    Parse MFE from RNAfold output.

    Compatible examples:
      ....((....)).... (-3.40)
      ....((....)).... ( -3.40)
      ....((....)).... ( -3.40 )
    """
    lines = rnafold_stdout.strip().splitlines()

    if len(lines) < 2:
        raise RuntimeError("Unexpected RNAfold output: " + rnafold_stdout)

    structure_line = lines[-1]

    match = re.search(r"\(\s*([-+]?\d+(?:\.\d+)?)\s*\)", structure_line)

    if not match:
        raise RuntimeError("Cannot parse MFE from RNAfold output: " + structure_line)

    return float(match.group(1))


def run_rnafold(seq, rnafold_exe, noPS_arg):
    """
    Run RNAfold and return MFE.
    DNA sequence will be converted to RNA sequence before folding.
    """
    rna_seq = seq.upper().replace("T", "U")

    result = subprocess.run(
        [rnafold_exe, noPS_arg],
        input=rna_seq + "\n",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip())

    mfe = parse_mfe_from_rnafold_output(result.stdout)

    return mfe


def process_one_record(record, rnafold_exe, noPS_arg):
    seq_id, seq = record

    seq = (
        seq.replace(" ", "")
           .replace("\t", "")
           .replace("\r", "")
           .replace("\n", "")
           .upper()
    )

    length = len(seq)

    if length == 0:
        return {
            "seq_id": seq_id,
            "length": 0,
            "GC_content": None,
            "NMFE": None
        }

    gc_content = calc_gc_content(seq)

    try:
        mfe = run_rnafold(seq, rnafold_exe, noPS_arg)
        nmfe = mfe / length
    except Exception as e:
        sys.stderr.write(f"Warning: failed to calculate MFE for {seq_id}: {e}\n")
        nmfe = None

    return {
        "seq_id": seq_id,
        "length": length,
        "GC_content": gc_content,
        "NMFE": nmfe
    }


def format_value(x, digits=6):
    if x is None:
        return "NA"

    if isinstance(x, str):
        return x

    return f"{x:.{digits}f}"


def main():
    args = parse_args()

    if args.threads < 1:
        sys.exit("Error: --threads must be >= 1")

    noPS_arg = detect_rnafold_noPS_arg(args.rnafold)

    sys.stderr.write(f"Using RNAfold: {args.rnafold}\n")
    sys.stderr.write(f"Using RNAfold noPS argument: {noPS_arg}\n")

    check_rnafold(args.rnafold, noPS_arg)

    records = list(read_fasta(args.fasta))

    if len(records) == 0:
        sys.exit("Error: no FASTA records found.")

    with open(args.out, "w") as out:
        out.write("seq_id\tlength\tGC_content\tNMFE\n")

        if args.threads == 1:
            for record in records:
                res = process_one_record(record, args.rnafold, noPS_arg)

                out.write(
                    f"{res['seq_id']}\t"
                    f"{res['length']}\t"
                    f"{format_value(res['GC_content'])}\t"
                    f"{format_value(res['NMFE'])}\n"
                )

        else:
            with ThreadPoolExecutor(max_workers=args.threads) as executor:
                results = executor.map(
                    lambda record: process_one_record(record, args.rnafold, noPS_arg),
                    records
                )

                for res in results:
                    out.write(
                        f"{res['seq_id']}\t"
                        f"{res['length']}\t"
                        f"{format_value(res['GC_content'])}\t"
                        f"{format_value(res['NMFE'])}\n"
                    )

    sys.stderr.write(f"Finished. Output written to: {args.out}\n")


if __name__ == "__main__":
    main()
