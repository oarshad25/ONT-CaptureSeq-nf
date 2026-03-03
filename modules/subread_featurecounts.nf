// take sorted and indexed bams and count reads overlapping genes
process SUBREAD_FEATURECOUNTS {
    label "low"

    publishDir "${params.outdir}/featurecounts/", mode: "copy"

    conda "bioconda subread=2.0.6"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/subread:2.0.6--he4a0461_2'
        : 'quay.io/biocontainers/subread:2.0.6--he4a0461_2'}"

    input:
    // sorted and indexed bam file
    tuple val(meta), path(bam), path(bai)

    // annotation gtf
    path annotation

    output:
    tuple val(meta), path("*geneCounts.txt"), emit: counts
    tuple val(meta), path("*geneCounts.txt.summary"), emit: summary

    script:
    """
    featureCounts \\
        -a ${annotation} \\
        -o "${meta.id}.geneCounts.txt" \\
        -L \\
        -t exon \\
        -g gene_id \\
        --primary \\
        --fraction \\
        -O \\
        -T ${task.cpus} \\
        ${bam}
    """
}
