// take sorted and indexed bams and count reads overlapping genes
process SUBREAD_FEATURECOUNTS {
    label "medium"

    publishDir "${params.outdir}/featurecounts/", mode: "copy"

    conda "bioconda subread=2.0.6"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/subread:2.0.6--he4a0461_2'
        : 'quay.io/biocontainers/subread:2.0.6--he4a0461_2'}"

    input:
    // list of sorted bam files
    path bams
    // list of corresponding bam indexes
    path bais
    // annotation gtf
    path annotation

    output:
    path "gene_counts.txt"

    script:
    """
    featureCounts \\
        -a ${annotation} \\
        -o gene_counts.txt \\
        -L \\
        -t exon \\
        -g gene_id \\
        --primary \\
        --fraction \\
        -O \\
        -T ${task.cpus} \\
        ${bams}
    """
}
