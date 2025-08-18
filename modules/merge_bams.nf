/*
* Use samtools to merge list of input BAMs and index merged bam
*/

process MERGE_BAMS {
    label 'low'

    conda "bioconda samtools=1.21"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.21--h96c455f_1'
        : 'quay.io/biocontainers/samtools:1.21--h96c455f_1'}"

    input:
    // list of bam files
    path bams
    // list of bam indexes
    path bais

    output:
    // merged and indexed bam file
    tuple path("merged.bam"), path("merged.bam.bai"), emit: merged_bambai

    script:
    """
    samtools merge --threads ${task.cpus} merged.bam ${bams}
    samtools index merged.bam
    """
}
