/*
* Use seqkit bam to generate tables of detailed BAM statistics
*/

process SEQKIT_BAM {
    tag "${meta.id}"
    label 'low'

    publishDir "${params.outdir}/qc/aligned/seqkit_bam/", mode: "copy"

    conda "bioconda samtools=2.9.0"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    input:
    // sorted bam file
    tuple val(meta), path(bam), path(bai)

    output:
    // table of detailed BAM statistics
    tuple val(meta), path("${meta.id}_seqkit_bam_stats.tsv"), emit: seqkit_bam_stats
    // histogram of aligned reference lengths
    tuple val(meta), path("${meta.id}_aln_len.pdf")
    // histogram of alignment accuracy
    tuple val(meta), path("${meta.id}_aln_acc.pdf")

    script:
    """
    # alignment length histogram
    seqkit bam -O ${meta.id}_aln_len.pdf -f RefAln ${bam}

    # alignment quality histogram
    seqkit bam -O ${meta.id}_aln_acc.pdf -f Acc ${bam}

    # Statistics from the sorted BAM:
    seqkit bam -s -Q ${bam} > ${meta.id}_seqkit_bam_stats.tsv 2>&1
    """
}
