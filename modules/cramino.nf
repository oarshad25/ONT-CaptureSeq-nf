/*
* Use cramino for quality assessment of alignment bams
*/

process CRAMINO {
    tag "${meta.id}"
    label 'low'

    publishDir "${params.outdir}/qc/aligned/cramino/", mode: "copy"

    conda "bioconda cramino=1.1.0"
    container "quay.io/biocontainers/cramino:1.1.0--h3dc2dae_0"

    input:
    // sorted bam file
    tuple val(meta), path(bam), path(bai)

    output:
    // BAM statistics
    tuple val(meta), path("${meta.id}_cramino_stats.txt"), emit: stats_txt
    // histograms
    tuple val(meta), path("${meta.id}_cramino_histograms.txt"), emit: histograms_txt

    script:
    """
    cramino \\
        --threads ${task.cpus} \\
        --hist ${meta.id}_cramino_histograms.txt \\
        --scaled \\
        --spliced \\
        --format 'text' \\
        ${bam} \\
        > ${meta.id}_cramino_stats.txt
    """
}
