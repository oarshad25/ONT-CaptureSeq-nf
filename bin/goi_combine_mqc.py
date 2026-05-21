#!/usr/bin/env python3
# bin/goi_combine_mqc.py

"""
goi_combine_mqc.py

Merges per-GOI metric TSVs from goi_readcount_n50.py, goi_fulllength_count.py,
and optionally goi_coverage_summary.py into a single MultiQC-compatible table.

Each input TSV is wide-format with one row per sample. This script merges
them horizontally by sample and prepends the MultiQC custom content header,
producing a table section in the MultiQC report.

Collision detection ensures that duplicate column names across input files
(which would indicate non-unique gene names in the GOI BED) cause an early
error rather than silently overwriting values.

The final MultiQC table is presented gene-first: columns are renamed to
<gene>__<metric> and sorted alphabetically by gene, with all available metrics
for a gene grouped together.

Usage:
    goi_combine_mqc.py --readcount-n50 sample_goi_readcount_n50.tsv
                       --fulllength sample_goi_fulllength.tsv
                       --coverage sample_goi_coverage.tsv
                       --sample SAMPLE_ID
                       --output sample_goi_summary_mqc.tsv
"""

import argparse
import csv
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine per-GOI metric TSVs into a single MultiQC table"
    )
    parser.add_argument("--readcount-n50", required=True, dest="readcount_n50",
                        help="Output of goi_readcount_n50.py")
    parser.add_argument("--fulllength",    required=True,
                        help="Output of goi_fulllength_count.py")
    parser.add_argument("--coverage",      required=False,
                        help="Optional output of goi_coverage_summary.py")
    parser.add_argument("--sample",        required=True,
                        help="Sample name")
    parser.add_argument("--output",        required=True,
                        help="Output MultiQC TSV path")
    return parser.parse_args()


# MultiQC custom content header template.
MQC_HEADER_TEMPLATE = """\
# id: 'goi_summary'
# section_name: 'Genes of Interest — Summary'
# description: 'Read count, percentage of total primary mapped reads, read N50, transcript-boundary-spanning read count, mean exonic depth, and exonic breadth per gene of interest locus.'
# format: 'tsv'
# plot_type: 'table'
# pconfig:
#   id: 'goi_summary_table'
#   title: 'Genes of Interest Summary'
"""


METRIC_PREFIXES = [
    ("read_count", "read_count_"),
    ("pct_total_primary_mapped", "pct_total_primary_mapped_"),
    ("n50", "n50_"),
    ("fulllength", "fulllength_"),
    ("mean_exonic_depth", "mean_exonic_depth_"),
    ("breadth_1x_pct", "breadth_1x_pct_"),
]


def write_column_headers(out, cols):
    """
    Write MultiQC table header config in the same order as the TSV columns.

    MultiQC tables default to one decimal place, which hides small
    percentages. Defining all headers preserves the input column order.
    """
    out.write("# headers:\n")
    for col in cols:
        metric = metric_from_output_column(col)
        out.write(f"#   {col}:\n")
        out.write(f"#     title: '{col}'\n")
        if metric == "pct_total_primary_mapped":
            out.write("#     format: '{:,.4f}'\n")
            out.write("#     suffix: '%'\n")
        elif metric == "breadth_1x_pct":
            out.write("#     format: '{:,.2f}'\n")
            out.write("#     suffix: '%'\n")
        elif metric == "mean_exonic_depth":
            out.write("#     format: '{:,.2f}'\n")
        else:
            out.write("#     format: '{:,.0f}'\n")


def split_input_metric_column(col):
    """
    Split an intermediate GOI metric column into metric and gene.

    Intermediate TSVs intentionally keep their existing metric-first schema
    (<metric>_<gene>). Gene names can contain underscores, so parse using the
    known metric prefixes rather than splitting on underscores.
    """
    for metric, prefix in METRIC_PREFIXES:
        if col.startswith(prefix):
            gene = col[len(prefix):]
            if gene:
                return metric, gene
    return None, None


def metric_from_output_column(col):
    """
    Return the metric suffix from a final gene-first column name.

    Final MultiQC columns use <gene>__<metric>. rsplit keeps this robust if a
    gene name itself contains the double-underscore delimiter.
    """
    if "__" in col:
        return col.rsplit("__", 1)[1]

    # Fallback for any future unrenamed recognised columns.
    metric, _gene = split_input_metric_column(col)
    return metric


def make_gene_first_table_row(merged):
    """
    Convert merged metric-first input columns to gene-first output columns.

    Returns an ordered dict-like mapping of final column name -> value. Known
    GOI metric columns are grouped by gene and metric order. Unknown columns are
    appended unchanged at the end so future metrics remain visible.
    """
    gene_values = {}
    unknown_values = {}

    for col, value in merged.items():
        if col == "Sample":
            continue

        metric, gene = split_input_metric_column(col)
        if metric is None:
            unknown_values[col] = value
            continue

        gene_values.setdefault(gene, {})[metric] = value

    final_values = {}
    for gene in sorted(gene_values, key=lambda g: (g.casefold(), g)):
        for metric, _prefix in METRIC_PREFIXES:
            if metric not in gene_values[gene]:
                continue
            final_col = f"{gene}__{metric}"
            if final_col in final_values or final_col in unknown_values:
                raise ValueError(f"Duplicate final output column '{final_col}'.")
            final_values[final_col] = gene_values[gene][metric]

    final_values.update(unknown_values)
    return final_values


def read_wide_tsv(path):
    """
    Read a wide-format TSV file and return the single data row as a dict.

    Each per-sample output TSV from the GOI metric scripts contains exactly
    one data row (one row per sample). Raises ValueError if this expectation
    is not met, catching truncated or malformed files early.

    Args:
        path: Path to the TSV file.

    Returns:
        Dict of column name -> value for the single data row.

    Raises:
        ValueError: If the file does not contain exactly one data row.
    """
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows   = list(reader)
    if len(rows) != 1:
        raise ValueError(
            f"Expected exactly 1 data row in {path}, got {len(rows)}."
        )
    return rows[0]


def validate_sample(row, path, expected_sample):
    """
    Warn if the sample name embedded in an input TSV does not match the
    expected sample name passed via --sample.

    A mismatch indicates that the wrong file has been passed for this sample,
    which could happen if Nextflow staging produces unexpected filename
    collisions.

    Args:
        row:             Dict of the data row from the TSV.
        path:            Path to the TSV (used in the warning message).
        expected_sample: Sample name from the --sample argument.
    """
    if row.get("Sample") != expected_sample:
        print(
            f"Warning: sample name in {path} ({row.get('Sample')}) "
            f"does not match expected sample ({expected_sample})",
            file=sys.stderr
        )


def check_collisions(named_dicts):
    """
    Check for duplicate column names across the input TSV dicts.

    Duplicate column names would arise if two GOIs share the same name in
    the GOI BED file, causing silent value overwriting during dict merge.
    Raises ValueError with a diagnostic message identifying which files
    share the conflicting column.

    Args:
        named_dicts: List of (name, dict) tuples where name is a label for
                     the input file and dict is the data row.

    Raises:
        ValueError: If any non-Sample column name appears in more than one input.
    """
    seen = {}
    for name, d in named_dicts:
        for col in d:
            if col == "Sample":
                continue
            if col in seen:
                raise ValueError(
                    f"Column '{col}' appears in both {seen[col]} and {name}. "
                    "Check that gene names are unique in your GOI BED file."
                )
            seen[col] = name


def main():
    args = parse_args()

    inputs = {
        "readcount_n50": read_wide_tsv(args.readcount_n50),
        "fulllength":    read_wide_tsv(args.fulllength),
    }
    if args.coverage:
        inputs["coverage"] = read_wide_tsv(args.coverage)

    # Validate sample names match across all input files
    for name, row in inputs.items():
        validate_sample(row, name, args.sample)

    # Check for column name collisions before merging
    check_collisions(list(inputs.items()))

    # Merge inputs — Sample column taken from readcount_n50, excluded from others
    merged = {**inputs["readcount_n50"]}
    for name in ["fulllength", "coverage"]:
        if name in inputs:
            merged.update({k: v for k, v in inputs[name].items() if k != "Sample"})

    final_values = make_gene_first_table_row(merged)
    cols = list(final_values)

    with open(args.output, "w") as out:
        out.write(MQC_HEADER_TEMPLATE)
        write_column_headers(out, cols)
        out.write("Sample\t"    + "\t".join(cols) + "\n")
        out.write(args.sample + "\t" + "\t".join(final_values[c] for c in cols) + "\n")


if __name__ == "__main__":
    main()
