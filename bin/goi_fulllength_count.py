#!/usr/bin/env python3
# bin/goi_fulllength_count.py

"""
goi_fulllength_count.py

Computes the count of reads spanning transcript boundaries per gene of interest (GOI) from an
aligned BAM file and a GTF annotation.

For each GOI, the longest annotated transcript is identified from the GTF.
Reads overlapping the GOI region are then classified as transcript-boundary-spanning if their
aligned reference span covers the transcript boundaries within a configurable
tolerance, rather than requiring an exact end-to-end match.

Reference span is determined using pysam's get_reference_positions(), which
correctly handles all CIGAR operators including hard clipping, = and X
operators, and intron-spanning gaps in spliced alignments. This is more
robust than manual CIGAR string parsing.

An optional --margin parameter allows reads to start or end within a
tolerance of the transcript boundaries, accommodating minor positional
variation at read ends common in Nanopore data due to adapter trimming or
homopolymer errors.

Ensembl version suffixes are stripped from GTF gene_id and gene_name
attributes before matching, allowing unversioned gene lists to match
versioned GTF annotations, consistent with make_goi_bed.py.

Outputs a single wide-format TSV with one row per sample and a
fulllength_<gene> count column per GOI, suitable for ingestion by
goi_combine_mqc.py.

Usage:
    goi_fulllength_count.py --bam sample.bam --goi goi.bed --gtf annotation.gtf
                           --sample SAMPLE_ID --output sample_goi_fulllength.tsv
"""

import argparse
import pysam
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute transcript-boundary-spanning read count per gene of interest"
    )
    parser.add_argument("--bam",    required=True, help="Indexed BAM file")
    parser.add_argument("--goi",    required=True, help="BED file of genes of interest (col4 = gene name)")
    parser.add_argument("--gtf",    required=True, help="GTF annotation file")
    parser.add_argument("--sample", required=True, help="Sample name for output header")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument(
        "--margin", type=int, default=100,
        help="Allow reads to start/end within this many bp of transcript boundaries (default: 100)"
    )
    args = parser.parse_args()
    if args.margin < 0:
        parser.error("--margin must be >= 0")
    return args


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


def parse_attribute(attr_string, key):
    """
    Extract a quoted attribute value from a GTF attribute column string.

    Args:
        attr_string: The full attribute column string (GTF field 9).
        key:         The attribute key to search for (e.g. 'gene_name').

    Returns:
        The unquoted value string if found, otherwise None.
    """
    match = re.search(rf'{key} "([^"]+)"', attr_string)
    return match.group(1) if match else None


def strip_ensembl_version(identifier):
    """
    Remove the version suffix from an Ensembl gene ID or name.

    Ensembl IDs in GTF files are often versioned (e.g. ENSG00000012048.5)
    while user-supplied gene lists typically contain unversioned IDs
    (e.g. ENSG00000012048). Stripping the suffix allows these to match.

    The strip is gated on the Ensembl ID pattern (^ENS[A-Z]*\\d+\\.\\d+$)
    so gene symbols or other identifiers that happen to contain a dot
    are left unchanged.

    Args:
        identifier: A gene ID or name string from a GTF attribute field.

    Returns:
        The identifier with version suffix removed if it matches the
        Ensembl pattern, otherwise the original identifier unchanged.
    """
    if re.match(r'^ENS[A-Z]*\d+\.\d+$', identifier):
        return identifier.rsplit('.', 1)[0]
    return identifier


def get_longest_transcripts(gtf_path, gene_names):
    """
    Identify the longest annotated transcript for each GOI gene from a GTF.

    Longest is defined as the greatest span from transcript start to end,
    i.e. the genomic footprint of the transcript including introns. This is
    used as the reference span that a transcript-boundary-spanning read must cover.

    Matching uses gene_name first, then gene_id as fallback, consistent with
    the approach in make_genes_bed.py. Ensembl version suffixes are stripped
    from both gene_id and gene_name GTF attributes before matching, allowing
    unversioned gene lists to match versioned GTF annotations.

    Args:
        gtf_path:   Path to the GTF annotation file.
        gene_names: Iterable of gene name/ID strings to match against.

    Returns:
        Dict mapping gene name/ID -> (chrom, tx_start, tx_end) where
        coordinates are 0-based half-open (BED-like).
    """
    longest  = {}
    gene_set = set(gene_names)

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")

            # Only process transcript-level features for boundary coordinates
            if len(fields) < 9 or fields[2] != "transcript":
                continue

            attr      = fields[8]
            gene_name = parse_attribute(attr, "gene_name")
            gene_id   = parse_attribute(attr, "gene_id")

            # Strip Ensembl version suffixes from GTF attributes before matching,
            # consistent with make_genes_bed.py — gene lists are assumed to be
            # unversioned while GTF annotations frequently carry version suffixes
            gene_name_stripped = strip_ensembl_version(gene_name) if gene_name else None
            gene_id_stripped   = strip_ensembl_version(gene_id)   if gene_id   else None

            # Match on stripped gene_name first, fall back to stripped gene_id
            matched = None
            if gene_name_stripped and gene_name_stripped in gene_set:
                matched = gene_name_stripped
            elif gene_id_stripped and gene_id_stripped in gene_set:
                matched = gene_id_stripped

            if matched is None:
                continue

            # Convert GTF 1-based inclusive coordinates to 0-based half-open
            chrom  = fields[0]
            start  = int(fields[3]) - 1
            end    = int(fields[4])
            length = end - start

            # Keep only the longest transcript per gene
            if matched not in longest or length > (longest[matched][2] - longest[matched][1]):
                longest[matched] = (chrom, start, end)

    return longest


def compute_fulllength_count(bam_path, chrom, region_start, region_end,
                             tx_start, tx_end, margin):
    """
    Compute the number of reads overlapping a GOI region whose aligned
    reference span covers the annotated start and end of its longest transcript.

    A read is classified as transcript-boundary-spanning if its aligned reference span starts
    at or before (tx_start + margin) and ends at or after (tx_end - margin).
    This relaxes the definition to tolerate modest 5' and 3' end truncation,
    which is common in Nanopore cDNA data.

    Reference span is derived from pysam's get_reference_positions(), which
    returns the reference coordinate of every aligned base, correctly handling
    all CIGAR operators. Reads with no aligned positions (e.g. fully
    soft-clipped) are counted in the total but not as transcript-boundary-spanning.

    Secondary and supplementary alignments are excluded to avoid
    double-counting chimeric or multi-mapped Nanopore reads.

    Args:
        bam_path:     Path to the indexed BAM file.
        chrom:        Chromosome/contig name.
        region_start: 0-based start of the GOI region to fetch reads from.
        region_end:   0-based end of the GOI region to fetch reads from.
        tx_start:     0-based start of the longest transcript.
        tx_end:       0-based end of the longest transcript.
        margin:       Tolerance in bp at transcript boundaries.

    Returns:
        Tuple of (full_length_count, total_count).
    """
    total      = 0
    fulllength = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, region_start, region_end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            total += 1

            ref_positions = read.get_reference_positions()
            if not ref_positions:
                # Read has no aligned bases — count in total but not as transcript-boundary-spanning
                continue

            read_ref_start = ref_positions[0]
            read_ref_end   = ref_positions[-1] + 1  # convert to half-open

            if read_ref_start <= (tx_start + margin) and read_ref_end >= (tx_end - margin):
                fulllength += 1

    return fulllength, total


def main():
    args    = parse_args()
    goi     = parse_goi_bed(args.goi)
    genes   = [g[3] for g in goi]
    longest = get_longest_transcripts(args.gtf, genes)

    # Warn for GOIs with no transcript annotation — these will be reported as NA
    missing = set(genes) - set(longest.keys())
    if missing:
        print(
            f"Warning: {len(missing)} GOI gene(s) not found in GTF as transcript "
            f"features: {', '.join(sorted(missing))}",
            file=sys.stderr
        )

    results = []
    for chrom, start, end, gene in goi:
        if gene not in longest:
            # Cannot compute the transcript-boundary-spanning count without transcript boundary coordinates
            results.append((gene, "NA"))
            continue

        tx_chrom, tx_start, tx_end = longest[gene]
        if tx_chrom != chrom:
            print(
                f"Warning: GOI BED region for {gene} is on {chrom} but longest "
                f"transcript annotation is on {tx_chrom}; reporting NA.",
                file=sys.stderr
            )
            results.append((gene, "NA"))
            continue

        fulllength, total = compute_fulllength_count(
            args.bam, chrom, start, end,
            tx_start, tx_end, args.margin
        )

        results.append((gene, fulllength))

    # Write wide-format TSV — one row per sample, one transcript-boundary-spanning count column per GOI
    with open(args.output, "w") as out:
        out.write(
            "Sample\t"
            + "\t".join(f"fulllength_{g}" for g, _ in results)
            + "\n"
        )
        out.write(
            args.sample + "\t"
            + "\t".join(str(r) for _, r in results)
            + "\n"
        )


if __name__ == "__main__":
    main()
