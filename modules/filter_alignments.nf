/*
* Use samtools to filter out from input BAM (alignments):
* - Secondary alignments
* - Supplementary alignments
* - Unmapped reads
*/

process FILTER_ALIGNMENTS {
    tag "${meta.id}"
    label 'low'

    publishDir "${params.outdir}/bam/filtered/", mode: "link"

    conda "bioconda samtools=1.21"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.21--h96c455f_1'
        : 'quay.io/biocontainers/samtools:1.21--h96c455f_1'}"

    input:
    // bam file
    tuple val(meta), path(bam), path(bai)

    output:
    // filtered and indexed bam file (input BAM with unmapped reads filtered out)
    tuple val(meta), path("${meta.id}.filtered.bam"), path("${meta.id}.filtered.bam.bai"), emit: filtered_bambai

    script:
    """
    # filter out non-primary (secondary) and supplementary alignments and unmapped reads
    samtools view -@ ${task.cpus} -bh -F 2308 ${bam} | \\
    samtools sort -@ ${task.cpus} -o ${meta.id}.filtered.bam

    # index filtered BAM
    samtools index -@ ${task.cpus} ${meta.id}.filtered.bam
    """
}
