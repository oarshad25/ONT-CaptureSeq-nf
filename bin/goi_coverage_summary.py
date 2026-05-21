#!/usr/bin/env python3
# bin/goi_coverage_summary.py

"""
goi_coverage_summary.py

Summarise mosdepth exon coverage outputs into one wide-format TSV per sample.

Inputs are the mosdepth region mean-depth file, the mosdepth threshold file
created with --thresholds 1, and the merged GOI exon BED used as mosdepth's
--by input. The exon BED is used as the source of truth for gene order and
total exonic bases.

Usage:
    goi_coverage_summary.py --regions sample.goi_exons.regions.bed.gz
                            --thresholds sample.goi_exons.thresholds.bed.gz
                            --exons goi_exons.bed
                            --sample SAMPLE_ID
                            --output sample_goi_coverage.tsv
"""

import argparse
import csv
import gzip
from collections import OrderedDict, defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarise exon-level mosdepth coverage by GOI gene"
    )
    parser.add_argument("--regions",    required=True,
                        help="mosdepth *.regions.bed.gz file")
    parser.add_argument("--thresholds", required=True,
                        help="mosdepth *.thresholds.bed.gz file from --thresholds 1")
    parser.add_argument("--exons",      required=True,
                        help="Merged exon BED used as mosdepth --by input")
    parser.add_argument("--sample",     required=True,
                        help="Sample name")
    parser.add_argument("--output",     required=True,
                        help="Output wide-format TSV")
    return parser.parse_args()


def open_text(path):
    """
    Open plain text or gzip-compressed files transparently.

    mosdepth writes BED-like outputs as .gz files, while the exon BED is plain
    text. A single helper keeps the parsing functions simple.
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_exon_bases(path):
    """
    Read the merged GOI exon BED and return total exonic bases per gene.

    The BED name column is expected to contain the gene name. Because the BED
    has already been merged per gene, summing interval lengths gives the
    denominator for both weighted mean depth and breadth.
    """
    exon_bases = OrderedDict()
    with open_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"Expected at least 4 BED columns in {path}: {line}")

            start = int(fields[1])
            end = int(fields[2])
            gene = fields[3]
            length = end - start
            if length <= 0:
                raise ValueError(
                    f"Found non-positive interval length in {path}: {line}"
                )

            exon_bases.setdefault(gene, 0)
            exon_bases[gene] += length

    if not exon_bases:
        raise ValueError(f"No exon intervals found in {path}.")

    return exon_bases


def read_region_depth_bases(path):
    """
    Read mosdepth region means and return sum(mean_depth * interval_length).

    mosdepth may preserve extra BED columns, so this parser takes the gene from
    column 4 and the mean depth from the final column.
    """
    depth_bases = defaultdict(float)
    with open_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(
                    f"Expected mosdepth region BED plus mean depth in {path}: {line}"
                )

            start = int(fields[1])
            end = int(fields[2])
            gene = fields[3]
            mean_depth = float(fields[-1])
            depth_bases[gene] += mean_depth * (end - start)

    return depth_bases


def read_threshold_bases(path):
    """
    Read mosdepth threshold counts and return bases covered at >=1x per gene.

    The process calls mosdepth with a single threshold, --thresholds 1, so the
    final column is the number of bases in the interval covered by at least one
    aligned primary read.
    """
    covered_bases = defaultdict(int)
    with open_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(
                    f"Expected mosdepth threshold BED plus count in {path}: {line}"
                )

            gene = fields[3]
            covered_bases[gene] += int(float(fields[-1]))

    return covered_bases


def main():
    args = parse_args()

    exon_bases = read_exon_bases(args.exons)
    depth_bases = read_region_depth_bases(args.regions)
    covered_bases = read_threshold_bases(args.thresholds)

    results = []
    for gene, total_bases in exon_bases.items():
        # Use the exon BED denominator so genes with no observed coverage still
        # report as 0 depth / 0 breadth rather than disappearing from the table.
        mean_depth = depth_bases[gene] / total_bases
        breadth_pct = covered_bases[gene] / total_bases * 100
        results.append((gene, mean_depth, breadth_pct))

    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        header = ["Sample"]
        row = [args.sample]

        for gene, mean_depth, breadth_pct in results:
            header.extend([
                f"mean_exonic_depth_{gene}",
                f"breadth_1x_pct_{gene}",
            ])
            row.extend([f"{mean_depth:.4f}", f"{breadth_pct:.4f}"])

        writer.writerow(header)
        writer.writerow(row)


if __name__ == "__main__":
    main()
