#!/usr/bin/env python3

import argparse
import os


def get_sample_name(log_file):
    """
    自动识别样品名：
    1. CK_1.Log.final.out -> CK_1
    2. Log.final.out -> 使用上一级目录名
    """
    basename = os.path.basename(log_file)

    if basename.endswith(".Log.final.out"):
        sample = basename.replace(".Log.final.out", "")
        if sample:
            return sample

    return os.path.basename(os.path.dirname(log_file))


def clean_value(value):
    """
    去掉百分号和空格
    """
    return value.strip().replace("%", "")


def to_int(value):
    try:
        return int(value)
    except Exception:
        return 0


def to_float(value):
    try:
        return float(value)
    except Exception:
        return 0.0


def parse_star_log(log_file):
    sample = get_sample_name(log_file)

    raw = {}

    with open(log_file, "r") as f:
        for line in f:
            if "|" not in line:
                continue

            key, value = line.rstrip("\n").split("|", 1)
            key = key.strip()
            value = clean_value(value)

            raw[key] = value

    total_reads = to_int(raw.get("Number of input reads", 0))

    unique_reads = to_int(raw.get("Uniquely mapped reads number", 0))
    unique_rate = to_float(raw.get("Uniquely mapped reads %", 0))

    multi_reads = to_int(raw.get("Number of reads mapped to multiple loci", 0))
    multi_rate = to_float(raw.get("% of reads mapped to multiple loci", 0))

    too_many_reads = to_int(raw.get("Number of reads mapped to too many loci", 0))
    too_many_rate = to_float(raw.get("% of reads mapped to too many loci", 0))

    unmapped_mismatch_reads = to_int(raw.get("Number of reads unmapped: too many mismatches", 0))
    unmapped_mismatch_rate = to_float(raw.get("% of reads unmapped: too many mismatches", 0))

    unmapped_short_reads = to_int(raw.get("Number of reads unmapped: too short", 0))
    unmapped_short_rate = to_float(raw.get("% of reads unmapped: too short", 0))

    unmapped_other_reads = to_int(raw.get("Number of reads unmapped: other", 0))
    unmapped_other_rate = to_float(raw.get("% of reads unmapped: other", 0))

    multiple_mapped_reads = multi_reads + too_many_reads
    multiple_mapped_rate = multi_rate + too_many_rate

    unmapped_reads = (
        unmapped_mismatch_reads
        + unmapped_short_reads
        + unmapped_other_reads
    )

    unmapped_rate = (
        unmapped_mismatch_rate
        + unmapped_short_rate
        + unmapped_other_rate
    )

    return {
        "Sample": sample,
        "Total_reads": str(total_reads),
        "Unique_mapped_reads": str(unique_reads),
        "Unique_mapped_rate": f"{unique_rate:.2f}",
        "Multiple_mapped_reads": str(multiple_mapped_reads),
        "Multiple_mapped_rate": f"{multiple_mapped_rate:.2f}",
        "Unmapped_reads": str(unmapped_reads),
        "Unmapped_rate": f"{unmapped_rate:.2f}",
    }


def main():
    parser = argparse.ArgumentParser(
        description="Extract STAR mapping statistics from Log.final.out files."
    )

    parser.add_argument(
        "log_files",
        nargs="+",
        help="Input STAR Log.final.out files"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="STAR_mapping.stat",
        help="Output file, default: STAR_mapping.stat"
    )

    args = parser.parse_args()

    columns = [
        "Sample",
        "Total_reads",
        "Unique_mapped_reads",
        "Unique_mapped_rate",
        "Multiple_mapped_reads",
        "Multiple_mapped_rate",
        "Unmapped_reads",
        "Unmapped_rate",
    ]

    with open(args.output, "w") as out:
        out.write("\t".join(columns) + "\n")

        for log_file in args.log_files:
            stat = parse_star_log(log_file)
            out.write("\t".join(stat[col] for col in columns) + "\n")

    print(f"STAR mapping statistics written to: {args.output}")


if __name__ == "__main__":
    main()