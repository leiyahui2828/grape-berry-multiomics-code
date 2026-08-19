#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed


def open_fastq(file):
    """
    自动识别 gzip 压缩 fastq
    """
    if file.endswith(".gz"):
        return gzip.open(file, "rt")
    else:
        return open(file, "r")


def get_sample_name(fq):
    """
    根据文件名自动生成样品名
    """
    name = os.path.basename(fq)

    for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.endswith(suffix):
            name = name.replace(suffix, "")

    # 去掉常见双端测序后缀
    for tag in ["_R1", "_R2", ".R1", ".R2", "_1", "_2", "_r1", "_r2"]:
        if name.endswith(tag):
            name = name[:-len(tag)]

    return name


def count_length_distribution(fq):
    """
    统计单个 fastq 文件的 reads 长度分布
    """
    length_count = Counter()
    total_reads = 0

    with open_fastq(fq) as f:
        while True:
            header = f.readline()
            if not header:
                break

            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline()

            if not qual:
                break

            read_len = len(seq)
            length_count[read_len] += 1
            total_reads += 1

    return length_count, total_reads


def process_one_fastq(fq):
    """
    单个进程处理一个 fastq 文件
    """
    sample = get_sample_name(fq)
    length_count, total_reads = count_length_distribution(fq)

    return {
        "sample": sample,
        "fq": fq,
        "length_count": dict(length_count),
        "total_reads": total_reads
    }


def write_single_sample_table(sample, length_count, outdir):
    """
    输出单个样品长度分布表
    """
    outfile = os.path.join(outdir, f"{sample}.length_distribution.tsv")

    total = sum(length_count.values())

    with open(outfile, "w") as out:
        out.write("Length\tCount\tPercent\n")

        for length in sorted(length_count):
            count = length_count[length]
            percent = count / total * 100 if total > 0 else 0
            out.write(f"{length}\t{count}\t{percent:.4f}\n")


def summarize_sample(sample, length_count):
    """
    生成单个样品统计信息
    """
    total_reads = sum(length_count.values())

    if total_reads == 0:
        return {
            "Sample": sample,
            "Total_reads": 0,
            "Min_length": "NA",
            "Max_length": "NA",
            "Mean_length": "NA",
            "Peak_length": "NA",
            "Peak_count": "NA",
            "Peak_percent": "NA"
        }

    lengths = sorted(length_count.keys())

    total_bases = sum(length * count for length, count in length_count.items())
    mean_len = total_bases / total_reads

    peak_length, peak_count = max(length_count.items(), key=lambda x: x[1])
    peak_percent = peak_count / total_reads * 100

    return {
        "Sample": sample,
        "Total_reads": total_reads,
        "Min_length": min(lengths),
        "Max_length": max(lengths),
        "Mean_length": f"{mean_len:.2f}",
        "Peak_length": peak_length,
        "Peak_count": peak_count,
        "Peak_percent": f"{peak_percent:.4f}"
    }


def write_summary_table(summary_list, outfile):
    """
    输出样品总体统计表
    """
    headers = [
        "Sample",
        "Total_reads",
        "Min_length",
        "Max_length",
        "Mean_length",
        "Peak_length",
        "Peak_count",
        "Peak_percent"
    ]

    with open(outfile, "w") as out:
        out.write("\t".join(headers) + "\n")

        for item in summary_list:
            out.write("\t".join(str(item[h]) for h in headers) + "\n")


def write_matrix_table(all_counts, outfile):
    """
    输出所有样品的长度分布 count 矩阵
    """
    all_lengths = sorted(
        set(length for sample in all_counts for length in all_counts[sample])
    )

    samples = list(all_counts.keys())

    with open(outfile, "w") as out:
        out.write("Length\t" + "\t".join(samples) + "\n")

        for length in all_lengths:
            values = [str(all_counts[sample].get(length, 0)) for sample in samples]
            out.write(str(length) + "\t" + "\t".join(values) + "\n")


def write_percent_matrix_table(all_counts, outfile):
    """
    输出所有样品的长度分布 percent 矩阵
    """
    all_lengths = sorted(
        set(length for sample in all_counts for length in all_counts[sample])
    )

    samples = list(all_counts.keys())

    sample_totals = {
        sample: sum(all_counts[sample].values())
        for sample in samples
    }

    with open(outfile, "w") as out:
        out.write("Length\t" + "\t".join(samples) + "\n")

        for length in all_lengths:
            values = []
            for sample in samples:
                total = sample_totals[sample]
                count = all_counts[sample].get(length, 0)
                percent = count / total * 100 if total > 0 else 0
                values.append(f"{percent:.4f}")

            out.write(str(length) + "\t" + "\t".join(values) + "\n")


def plot_length_distribution(all_counts, outfile):
    """
    默认绘制所有样品的 reads 长度百分比分布图
    坐标轴刻度使用整数
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator
    except ImportError:
        print("Warning: matplotlib not installed, skip plotting.")
        return

    plt.figure(figsize=(8, 5))

    all_lengths = []

    for sample, length_count in all_counts.items():
        lengths = sorted(length_count.keys())
        counts = [length_count[x] for x in lengths]
        total = sum(counts)

        if total == 0:
            continue

        percents = [x / total * 100 for x in counts]
        all_lengths.extend(lengths)

        plt.plot(
            lengths,
            percents,
            marker="o",
            linewidth=1,
            markersize=3,
            label=sample
        )

    plt.xlabel("Read length (nt)")
    plt.ylabel("Reads percentage (%)")
    plt.title("FASTQ read length distribution")

    ax = plt.gca()

    # x 轴为 reads 长度，必须是整数
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    # y 轴百分比也显示为整数刻度
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # 可选：让 x 轴范围刚好覆盖 reads 长度
    if all_lengths:
        ax.set_xlim(min(all_lengths), max(all_lengths))

    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Calculate read length distribution from fq.gz / fastq.gz files."
    )

    parser.add_argument(
        "-i", "--input",
        nargs="+",
        required=True,
        help="Input fq.gz / fastq.gz files"
    )

    parser.add_argument(
        "-o", "--outdir",
        default="length_distribution",
        help="Output directory [default: length_distribution]"
    )

    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=1,
        help="Number of parallel processes [default: 1]"
    )

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    all_counts = {}
    summary_list = []

    max_workers = max(1, args.threads)

    print(f"Input FASTQ files: {len(args.input)}")
    print(f"Parallel processes: {max_workers}")

    results = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_fq = {
            executor.submit(process_one_fastq, fq): fq
            for fq in args.input
        }

        for future in as_completed(future_to_fq):
            fq = future_to_fq[future]

            try:
                result = future.result()
                results.append(result)
                print(f"Finished: {result['sample']}  {fq}")
            except Exception as e:
                print(f"Error processing {fq}: {e}")

    # 保持输出样品顺序与输入文件顺序一致
    sample_order = [get_sample_name(fq) for fq in args.input]
    result_dict = {x["sample"]: x for x in results}

    for sample in sample_order:
        if sample not in result_dict:
            continue

        result = result_dict[sample]
        length_count = result["length_count"]

        all_counts[sample] = length_count

        write_single_sample_table(
            sample,
            length_count,
            args.outdir
        )

        summary = summarize_sample(sample, length_count)
        summary_list.append(summary)

    write_summary_table(
        summary_list,
        os.path.join(args.outdir, "length_distribution.summary.tsv")
    )

    write_matrix_table(
        all_counts,
        os.path.join(args.outdir, "length_distribution.count_matrix.tsv")
    )

    write_percent_matrix_table(
        all_counts,
        os.path.join(args.outdir, "length_distribution.percent_matrix.tsv")
    )

    plot_length_distribution(
        all_counts,
        os.path.join(args.outdir, "length_distribution.pdf")
    )

    print("Done.")


if __name__ == "__main__":
    main()
