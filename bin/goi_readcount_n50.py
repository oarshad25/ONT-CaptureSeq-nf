#!/usr/bin/env python3
# bin/goi_readcount_n50.py

"""
goi_readcount_n50.py

Computes read count, percentage of total primary mapped reads, and read
length N50 for each gene of interest (GOI) region from an aligned BAM file.

For each region defined in the GOI BED file, primary mapped reads are
extracted and their query lengths collected. Query length (read.query_length)
includes soft-clipped bases but excludes hard-clipped bases, reflecting the
length of the sequence present in the BAM record.

Secondary and supplementary alignments are excluded to avoid double-counting
chimeric or multi-mapped Nanopore reads.

Outputs a single wide-format TSV with one row per sample and grouped
read_count_<gene>, pct_total_primary_mapped_<gene>, and n50_<gene> columns
for each GOI, suitable for ingestion by goi_combine_mqc.py.

Usage:
    goi_readcount_n50.py --bam sample.bam --goi goi.bed
                         --sample SAMPLE_ID --output sample_goi_readcount_n50.tsv
"""

import argparse
import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Count primary mapped reads, compute percentage of total primary mapped reads, and compute N50 of read lengths per GOI region"
    )
    parser.add_argument("--bam",    required=True, help="Indexed BAM file")
    parser.add_argument("--goi",    required=True, help="BED file of genes of interest (col4 = gene name)")
    parser.add_argument("--sample", required=True, help="Sample name for output header")
    parser.add_argument("--output", required=True, help="Output TSV path")
    return parser.parse_args()


def parse_goi_bed(path):
    """
    Parse a BED file of genes of interest into a list of region tuples.

    Skips blank lines and comment lines (starting with '#').
    Expects at least 4 columns: chrom, start (0-based), end, name.

    Args:
        path: Path to the GOI BED file.

    Returns:
        List of (chrom, start, end, gene_name) tuples.
    """
    regions = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            regions.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
    return regions


def compute_n50(lengths):
    """
    Compute the N50 read length from a list of read lengths.

    N50 is the length L such that reads of length >= L account for at
    least half of the total sequenced bases. It is computed by sorting
    lengths in descending order and finding the length at which the
    cumulative sum first reaches or exceeds half the total.

    For Nanopore CaptureSeq, a high N50 relative to the expected
    transcript length suggests reads may be long enough to support
    isoform-level interpretation.

    Args:
        lengths: List of integer read lengths.

    Returns:
        Integer N50 value, or 0 if the list is empty.
    """
    if not lengths:
        return 0
    sorted_lengths = sorted(lengths, reverse=True)
    total          = sum(sorted_lengths)
    cumsum         = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= total / 2:
            return length
    return 0


def get_readcount_and_n50(bam_path, chrom, start, end):
    """
    Extract read lengths for all primary mapped reads overlapping a region
    and compute the read count and N50.

    Excludes secondary alignments (flag 0x100) and supplementary alignments
    (flag 0x800) to avoid double-counting reads that have split or chimeric
    alignments, which are common in Nanopore long-read data.

    Args:
        bam_path: Path to the indexed BAM file.
        chrom:    Chromosome/contig name.
        start:    0-based region start coordinate.
        end:      0-based region end coordinate.

    Returns:
        Tuple of (read_count, n50) where read_count is the number of
        primary mapped reads and n50 is the N50 of their query lengths.
    """
    lengths = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            lengths.append(read.query_length)
    return len(lengths), compute_n50(lengths)


def get_total_primary_mapped_reads(bam_path):
    """
    Count all primary mapped reads in the BAM.

    Returns:
        Integer count of mapped reads excluding secondary and supplementary
        alignments.
    """
    count = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            count += 1
    return count


def main():
    args = parse_args()
    goi  = parse_goi_bed(args.goi)
    total_primary_mapped = get_total_primary_mapped_reads(args.bam)

    results = []
    for chrom, start, end, gene in goi:
        count, n50 = get_readcount_and_n50(args.bam, chrom, start, end)
        pct = (count / total_primary_mapped * 100) if total_primary_mapped else 0
        results.append((gene, count, pct, n50))

    # Write wide-format TSV — one row per sample, grouped columns per GOI:
    # read_count_<gene>, pct_total_primary_mapped_<gene>, n50_<gene>
    with open(args.output, "w") as out:
        out.write(
            "Sample\t"
            + "\t".join(
                f"read_count_{g}\tpct_total_primary_mapped_{g}\tn50_{g}"
                for g, _, _, _ in results
            )
            + "\n"
        )
        out.write(
            args.sample + "\t"
            + "\t".join(f"{c}\t{p:.4f}\t{n}" for _, c, p, n in results)
            + "\n"
        )


if __name__ == "__main__":
    main()
