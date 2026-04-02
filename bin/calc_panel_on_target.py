#!/usr/bin/env python
"""Calculate per-sample panel on-target metrics from a gene counts matrix.

This script reads a gene-by-sample count table, matches the gene identifier
column against a plain-text list of capture-panel targets, and writes a
MultiQC-compatible table containing total assigned reads, on-target reads, and
the percentage of reads assigned to panel genes for each sample.
"""

import argparse
import sys

import pandas as pd


def calculate_on_target(counts_file, panel_targets, output_file):
    """Summarise on-target counts for each sample in a counts matrix.

    Args:
        counts_file: Path to the tab-delimited gene-by-sample counts matrix.
        panel_targets: Path to a file containing one target gene ID per line.
        output_file: Path to the MultiQC-formatted output table to create.

    Raises:
        SystemExit: If the counts file cannot be parsed.
        ValueError: If none of the requested panel genes are present in the
            counts matrix gene column.
    """
    # Read the counts matrix while ignoring any metadata header
    # lines that may be prefixed with '#'.
    try:
        df = pd.read_csv(counts_file, sep='\t', comment='#')
    except Exception as e:
        print(f"Error reading counts file: {e}")
        sys.exit(1)

    # Treat the first column as the gene identifier column and all remaining
    # columns as per-sample count vectors.
    gene_col = df.columns[0]
    sample_cols = df.columns[1:]

    # Strip Ensembl-style version suffixes so the counts matrix and panel list
    # can be matched whether they contain versioned or unversioned IDs.
    df[gene_col] = df[gene_col].astype(str).apply(lambda x: x.split('.')[0])

    with open(panel_targets, 'r') as f:
        target_genes = set(line.strip().split('.')[0] for line in f if line.strip())

    # Precompute which count-matrix rows belong to the target panel so we can
    # reuse the match set for each sample and fail fast if the panel is absent.
    panel_gene_mask = df[gene_col].isin(target_genes)
    matched_targets = panel_gene_mask.any()
    if not matched_targets:
        raise ValueError(
            "No gene IDs from the panel target list were found in the "
            f"'{gene_col}' column of {counts_file}."
        )

    # Summarise absolute and percentage on-target counts for each sample using
    # the same precomputed panel membership mask.
    metrics = {}
    for sample in sample_cols:
        total_assigned = df[sample].sum()
        on_target_reads = df.loc[panel_gene_mask, sample].sum()
        on_target_rate = (on_target_reads / total_assigned * 100) if total_assigned > 0 else 0

        metrics[sample] = {
            'Total_Assigned_Reads': round(total_assigned, 2),
            'On_Target_Reads': round(on_target_reads, 2),
            'On_Target_Rate_Pct': round(on_target_rate, 2)
        }

    # Use the sample name as the first column so MultiQC renders one row per sample.
    results_df = pd.DataFrame.from_dict(metrics, orient='index')
    results_df.index.name = 'Sample'

    with open(output_file, 'w') as f:
        # Emit the MultiQC custom-content metadata block before the tabular data.
        f.write("# id: 'capture_on_target'\n")
        f.write("# section_name: 'Panel On-Target Rate'\n")
        f.write("# description: 'Reads mapping to genes in the target panel.'\n")
        f.write("# plot_type: 'table'\n")

        # --- MULTIQC COLUMN CONFIGURATION ---
        f.write("# headers:\n")
        f.write("#   Total_Assigned_Reads:\n")
        f.write("#     title: 'Total Assigned'\n")
        f.write("#     format: '{:,.0f}'\n")       # Add commas, strip all decimals

        f.write("#   On_Target_Reads:\n")
        f.write("#     title: 'On Target'\n")
        f.write("#     format: '{:,.0f}'\n")       # Add commas, strip all decimals

        f.write("#   On_Target_Rate_Pct:\n")
        f.write("#     title: 'On Target Rate'\n")
        f.write("#     format: '{:,.1f}'\n")       # Keep 1 decimal place
        f.write("#     suffix: '%'\n")             # Add the percentage symbol
        f.write("#     min: 0\n")                  # Set scale minimum
        f.write("#     max: 100\n")                # Set scale maximum
        # ------------------------------------
        results_df.to_csv(f, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Calculate per-sample panel on-target metrics from a tab-delimited "
            "gene counts matrix and write a MultiQC-compatible table."
        )
    )
    parser.add_argument(
        "-c",
        "--counts",
        required=True,
        help="Tab-delimited gene-by-sample counts matrix.",
    )
    parser.add_argument(
        "-t",
        "--targets",
        required=True,
        help="Text file listing panel gene IDs, one per line.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="on_target_mqc.tsv",
        help="Output path for the MultiQC table.",
    )
    args = parser.parse_args()

    calculate_on_target(args.counts, args.targets, args.output)
