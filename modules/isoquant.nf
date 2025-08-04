process ISOQUANT {
    label 'medium'

    publishDir "${params.outdir}", mode: 'link'

    conda "bioconda::isoquant=3.7.0"
    container "quay.io/biocontainers/isoquant:3.7.0--hdfd78af_0"

    input:
    // list of bam files
    path bams
    // list of bam indexes
    path bais
    // genome fasta file
    path genome
    // annotation gtf
    path annotation
    // Boolean, whether to set --complete_genedb option
    val complete_genedb

    output:
    path ("**/*.read_assignments.tsv.gz"), emit: read_assignments
    path ("**/*.corrected_reads.bed.gz"), emit: corrected_reads
    path ("**/*.transcript_counts.tsv"), emit: transcript_counts
    path ("**/*.gene_counts.tsv"), emit: gene_counts
    path ("**/*.transcript_tpm.tsv"), emit: transcript_tpm
    path ("**/*.gene_tpm.tsv"), emit: gene_tpm
    path ("**/isoquant.log"), emit: log
    path ("**/*.novel_vs_known.SQANTI-like.tsv"), emit: sqanti_output, optional: true
    path ("**/*.exon_counts.tsv"), emit: exon_counts, optional: true
    path ("**/*.intron_counts.tsv"), emit: intron_counts, optional: true
    path ("**/*.gene_grouped_counts.linear.tsv"), emit: gene_grouped_counts_linear
    path ("**/*.transcript_grouped_counts.linear.tsv"), emit: transcript_grouped_counts_linear
    path ("**/*.exon_grouped_counts.tsv"), emit: exon_grouped_counts, optional: true
    path ("**/*.intron_grouped_counts.tsv"), emit: intron_grouped_counts, optional: true
    path ("**/*.gene_grouped_counts.tsv"), emit: gene_grouped_counts, optional: true
    path ("**/*.transcript_grouped_counts.tsv"), emit: transcript_grouped_counts, optional: true
    path ("**/*.gene_grouped_tpm.tsv"), emit: gene_grouped_tpm, optional: true
    path ("**/*.transcript_grouped_tpm.tsv"), emit: transcript_grouped_tpm, optional: true
    path ("**/*.transcript_models.gtf"), emit: transcript_models, optional: true
    path ("**/*.transcript_model_reads.tsv.gz"), emit: transcript_model_reads, optional: true
    path ("**/*.discovered_transcript_counts.tsv"), emit: discovered_transcript_counts, optional: true
    path ("**/*.discovered_gene_counts.tsv"), emit: discovered_gene_counts, optional: true
    path ("**/*.discovered_transcript_tpm.tsv"), emit: discovered_transcript_tpm, optional: true
    path ("**/*.discovered_gene_tpm.tsv"), emit: discovered_gene_tpm, optional: true
    path ("**/*.extended_annotation.gtf"), emit: extended_annotation, optional: true
    path ("**/*.discovered_transcript_grouped_counts.linear.tsv"), emit: discovered_transcript_grouped_counts_linear, optional: true
    path ("**/*.discovered_gene_grouped_counts.linear.tsv"), emit: discovered_gene_grouped_counts_linear, optional: true
    path ("**/*.discovered_transcript_grouped_counts.tsv"), emit: discovered_transcript_grouped_counts, optional: true
    path ("**/*.discovered_gene_grouped_counts.tsv"), emit: discovered_gene_grouped_counts, optional: true

    script:
    def complete_genedb_flag = complete_genedb ? "--complete_genedb" : ""

    // use bam simpleName (meta.id) for sample labels
    def labels = bams.collect { fn -> "${fn.simpleName}" }.join(' ')

    """
    isoquant.py \\
        --threads ${task.cpus} \\
        --data_type nanopore \\
        --bam ${bams} \\
        --reference ${genome} \\
        --genedb ${annotation} ${complete_genedb_flag} \\
        --output isoquant \\
        --labels ${labels} \\
        --counts_format matrix
    """
}
