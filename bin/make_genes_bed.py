#!/usr/bin/env python3

"""
make_genes_bed.py

Creates a BED file of genomic regions for a list of genes of interest by
extracting coordinates from a GTF annotation file.

Matching strategy:
  - GTF gene_id field is version-stripped (e.g. ENSG00000012048.5 ->
    ENSG00000012048) before comparison, allowing gene lists containing
    unversioned Ensembl IDs to match versioned GTF annotations.
  - gene_name is tried as a fallback, allowing gene lists containing
    gene symbols (e.g. BRCA1) to match GTFs that carry both attributes.
  - Gene list entries are assumed to be clean (no version suffixes).

Outputs:
  - A sorted BED file of matched regions (--output).
  - If --merge-overlaps is set, overlapping or book-ended intervals are
    merged per gene, chromosome, and strand before writing the BED file.
  - A plain text file listing any genes from the input list that were
    not found in the GTF (--not-found). Only written if there are missing
    genes; absent if all genes were matched.

Usage:
    make_genes_bed.py --gtf annotation.gtf --gene-list genes.txt
                      --output genes.bed --not-found not_found.txt
"""

import argparse
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Create a BED file of GTF feature intervals for genes of interest "
            "from a GTF annotation and a plain text list of gene IDs or names. "
            "Use --feature exon --merge-overlaps to create merged exon BEDs "
            "for coverage calculations."
        )
    )
    parser.add_argument("--gtf",       required=True, help="Input GTF annotation file")
    parser.add_argument("--gene-list", required=True, help="Plain text file of gene IDs or gene names, one per line")
    parser.add_argument("--output",    required=True, help="Output BED file of matched feature intervals")
    parser.add_argument("--not-found", required=True, help="Output plain text file for requested genes with no matching feature in the GTF")
    parser.add_argument(
        "--feature", default="gene",
        help="GTF feature type to extract coordinates from (default: gene). "
             "Use 'exon' for exon intervals or 'transcript' for per-transcript intervals."
    )
    parser.add_argument(
        "--merge-overlaps", action="store_true",
        help="Merge overlapping or book-ended matched intervals per gene, "
             "chromosome, and strand before writing the BED file. Recommended "
             "for exon BEDs so exons shared by transcript isoforms are counted once."
    )
    return parser.parse_args()


def load_gene_list(path):
    """
    Load a plain text gene list into a set, ignoring blank lines and
    comments (lines starting with '#'). Duplicates are silently discarded
    by virtue of set membership.

    Args:
        path: Path to the gene list file.

    Returns:
        Set of gene ID/name strings.
    """
    genes = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            genes.add(line)
    return genes


def strip_ensembl_version(identifier):
    """
    Remove the version suffix from an Ensembl gene ID.

    Ensembl IDs in GTF files are often versioned (e.g. ENSG00000012048.5)
    while user-supplied gene lists typically contain unversioned IDs
    (e.g. ENSG00000012048). Stripping the suffix allows these to match.

    The strip is gated on the Ensembl ID pattern (^ENS[A-Z]*\\d+\\.\\d+$)
    so gene symbols or other identifiers that happen to contain a dot
    (e.g. a hypothetical symbol like GENE.2) are left unchanged.

    Args:
        identifier: A gene ID or name string from a GTF attribute field.

    Returns:
        The identifier with version suffix removed if it matches the
        Ensembl pattern, otherwise the original identifier unchanged.
    """
    if re.match(r'^ENS[A-Z]*\d+\.\d+$', identifier):
        return identifier.rsplit('.', 1)[0]
    return identifier


def parse_attribute(attr_string, key):
    """
    Extract a quoted attribute value from a GTF attribute column string.

    GTF attribute strings have the format:
        key1 "value1"; key2 "value2"; ...

    Args:
        attr_string: The full attribute column string (GTF field 9).
        key:         The attribute key to search for (e.g. 'gene_id').

    Returns:
        The unquoted value string if found, otherwise None.
    """
    match = re.search(rf'{key} "([^"]+)"', attr_string)
    return match.group(1) if match else None


def merge_overlapping_records(records):
    """
    Merge overlapping or book-ended BED records per gene/chromosome/strand.

    GTF exon features are usually repeated across transcript isoforms. Without
    this merge step, shared exon sequence would be counted more than once when
    calculating exon-level coverage breadth and depth.

    Args:
        records: BED records as (chrom, start, end, name, strand, matched).

    Returns:
        List of merged records in the same tuple format.
    """
    grouped = {}
    for chrom, start, end, name, strand, matched in records:
        # Keep genes separate even if their exons overlap on the genome; the
        # coverage summary is intended to be gene-level, not locus-level.
        key = (chrom, name, strand, matched)
        grouped.setdefault(key, []).append((start, end))

    merged_records = []
    for (chrom, name, strand, matched), intervals in grouped.items():
        intervals.sort()
        current_start, current_end = intervals[0]

        for start, end in intervals[1:]:
            # BED intervals are half-open, so start == current_end means the
            # intervals touch with no gap and can be treated as one region.
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                merged_records.append(
                    (chrom, current_start, current_end, name, strand, matched)
                )
                current_start, current_end = start, end

        merged_records.append(
            (chrom, current_start, current_end, name, strand, matched)
        )

    return merged_records


def main():
    args     = parse_args()
    gene_set = load_gene_list(args.gene_list)

    if not gene_set:
        sys.exit("Error: gene list is empty after filtering blank lines and comments.")

    found   = set()   # gene_set entries successfully matched in GTF
    records = []      # BED records to write: (chrom, start, end, name, strand, matched)

    with open(args.gtf) as f:
        for line in f:
            # Skip GTF header comment lines
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")

            # Skip malformed lines and lines for features other than the target type
            if len(fields) < 9:
                continue
            if fields[2] != args.feature:
                continue

            attr      = fields[8]
            gene_id   = parse_attribute(attr, "gene_id")
            gene_name = parse_attribute(attr, "gene_name")

            # Match against gene_set in priority order:
            #   1. Version-stripped gene_id  (handles versioned Ensembl GTFs)
            #   2. gene_name as-is           (handles gene symbol input lists)
            matched = None
            if gene_id and strip_ensembl_version(gene_id) in gene_set:
                matched = strip_ensembl_version(gene_id)
            elif gene_name and gene_name in gene_set:
                matched = gene_name

            if matched is None:
                continue

            # GTF coordinates are 1-based inclusive; convert start to 0-based for BED
            chrom  = fields[0]
            start  = int(fields[3]) - 1
            end    = int(fields[4])
            strand = fields[6]

            # Prefer gene_name as the BED name field for human readability;
            # fall back to gene_id if gene_name is absent
            name = gene_name if gene_name else gene_id

            records.append((chrom, start, end, name, strand, matched))
            found.add(matched)

    # Genes in the input list that had no matching feature in the GTF
    missing = gene_set - found

    # Only write the not-found file if there are missing genes — the calling
    # process declares it as an optional output and expects its absence when
    # all genes were successfully matched
    if missing:
        with open(args.not_found, "w") as out:
            for gene in sorted(missing):
                out.write(f"{gene}\n")
        print(
            f"Warning: {len(missing)} gene(s) not found in GTF as "
            f"'{args.feature}' features: {', '.join(sorted(missing))}",
            file=sys.stderr
        )

    if not records:
        sys.exit(
            "Error: no records written — no genes from the gene list were found "
            f"in the GTF as '{args.feature}' features. Check that gene IDs/names "
            "and the --feature type match your GTF."
        )

    if args.merge_overlaps:
        before_merge = len(records)
        records = merge_overlapping_records(records)
        print(
            f"Merged {before_merge} '{args.feature}' record(s) into "
            f"{len(records)} non-overlapping BED interval(s).",
            file=sys.stderr
        )

    # Sort by chrom then start position for consistent, predictable output
    records.sort(key=lambda r: (r[0], r[1]))

    with open(args.output, "w") as out:
        for chrom, start, end, name, strand, _ in records:
            out.write(f"{chrom}\t{start}\t{end}\t{name}\t.\t{strand}\n")

    print(
        f"Written {len(records)} record(s) for {len(found)} gene(s) to {args.output}. "
        f"{len(missing)} gene(s) not found — see {args.not_found}."
        if missing else
        f"Written {len(records)} record(s) for {len(found)} gene(s) to {args.output}. "
        f"All genes matched.",
        file=sys.stderr
    )


if __name__ == "__main__":
    main()
