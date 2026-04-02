#!/usr/bin/env python
"""Report per-sample read count and N50 metrics for select target genes.

This script reads one or more BAM files together with a BED file describing
gene intervals of interest, counts primary alignments overlapping each interval,
computes an N50 read length for those reads, and writes the results in a
MultiQC-compatible tabular format.
"""

import argparse
import os
import sys

import pysam


def calculate_n50(lengths):
    """Return the N50 value for a collection of read lengths."""
    if not lengths:
        return 0
    # Walk the lengths from longest to shortest until at least half of the
    # total sequenced bases are accounted for.
    lengths.sort(reverse=True)
    total_length = sum(lengths)
    half_total = total_length / 2.0
    running_sum = 0
    for length in lengths:
        running_sum += length
        if running_sum >= half_total:
            return length
    return 0

def calculate_read_stats(bam_files, bed_file, output_file):
    """Calculate per-gene read count and N50 metrics for each BAM file.

    Args:
        bam_files: Iterable of indexed BAM file paths to summarise.
        bed_file: BED file containing target intervals and output labels.
        output_file: Path to the MultiQC-formatted output table to create.
    """
    # Load the select-gene intervals and the label that should appear in MultiQC.
    bed_regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                bed_regions.append({
                    'chrom': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'label': parts[3]
                })

    # Keep output columns in the same order as the BED file and store metrics
    # in a nested sample -> gene-label mapping for easy row construction later.
    stats = {}
    labels = [r['label'] for r in bed_regions]
    # Track samples we had to skip so the task logs show incomplete inputs.
    dropped_samples = []

    for bam_path in bam_files:
        # Derive the sample ID from pipeline BAM names while preserving dots
        # that are part of the sample ID itself.
        sample_name = os.path.basename(bam_path)
        if sample_name.endswith(".bam"):
            sample_name = sample_name[:-4]
        if sample_name.endswith(".filtered"):
            sample_name = sample_name[:-9]
        stats[sample_name] = {}

        try:
            bam = pysam.AlignmentFile(bam_path, "rb")
        except Exception as e:
            # Record BAM-open failures and continue so other samples can still
            # be summarised in the same task.
            dropped_samples.append((sample_name, bam_path, str(e)))
            continue

        for row in bed_regions:
            read_lengths = []
            # Count primary alignments overlapping the gene interval and collect
            # their query lengths for the N50 calculation.
            for read in bam.fetch(row['chrom'], row['start'], row['end']):
                if not read.is_secondary and not read.is_supplementary:
                    read_lengths.append(read.query_length)

            # Record both depth and length-summary metrics for this sample/gene
            # pair so MultiQC can display them as adjacent columns.
            stats[sample_name][row['label']] = {
                'Read_Count': len(read_lengths),
                'N50_Length_bp': int(calculate_n50(read_lengths))
            }
        bam.close()

    if dropped_samples:
        # Emit skipped-sample details to stderr so they appear in the process log.
        for sample_name, bam_path, error_message in dropped_samples:
            print(
                f"Dropped sample '{sample_name}' from select-gene stats because "
                f"'{bam_path}' could not be opened: {error_message}",
                file=sys.stderr,
            )

    if not stats:
        print("No data collected.")
        return

    # Manually format the pivot table for MultiQC
    with open(output_file, 'w') as f:
        # Emit MultiQC custom-content metadata, then a wide table with one row
        # per sample and paired columns for each target gene.
        f.write("# id: 'specific_gene_stats'\n")
        f.write("# section_name: 'Specific Target Gene Statistics'\n")
        f.write("# description: 'Read counts and N50 for select genes.'\n")
        f.write("# plot_type: 'table'\n")

        # Write Header
        header = ["Sample"]
        for label in labels:
            header.extend([f"{label}_Read_Count", f"{label}_N50_Length_bp"])
        f.write('\t'.join(header) + '\n')

        # Write Rows
        for sample, gene_data in stats.items():
            row = [sample]
            for label in labels:
                row.extend([str(gene_data[label]['Read_Count']), str(gene_data[label]['N50_Length_bp'])])
            f.write('\t'.join(row) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Calculate per-sample read count and N50 statistics for selected "
            "genes from one or more indexed BAM files."
        )
    )
    parser.add_argument(
        "-b",
        "--bams",
        nargs='+',
        required=True,
        help="One or more indexed BAM files to summarise.",
    )
    parser.add_argument(
        "-bed",
        "--bed",
        required=True,
        help="BED file containing target intervals (chr, start, end as first three columns) and display labels (4th column).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="panel_metrics_gene_stats_mqc.tsv",
        help="Output path for the MultiQC table.",
    )
    args = parser.parse_args()
    calculate_read_stats(args.bams, args.bed, args.output)
