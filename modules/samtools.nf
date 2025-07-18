/*
* Use samtools to convert input sam to sorted bam, index it and
* compute alignment statistics with flagstat
*/

process SAMTOOLS {
    tag "${meta.id}"
    label 'low'

    publishDir "${params.outdir}/bam/", mode: "link"

    conda "bioconda samtools=1.21"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.21--h96c455f_1'
        : 'quay.io/biocontainers/samtools:1.21--h96c455f_1'}"

    input:
    // sam file
    tuple val(meta), path(sam)

    output:
    // sorted and indexed bam
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bai"), emit: bambai

    // flagstats
    tuple val(meta), path("${meta.id}.flagstat.txt"), emit: flagstat

    script:
    """
    # convert SAM to sorted BAM
    samtools view -Sb ${sam} | \\
    samtools sort -@ ${task.cpus} -O bam -o ${meta.id}.bam -

    # index BAM
    samtools index ${meta.id}.bam ${meta.id}.bai

    # post-mapping QC
    samtools flagstat ${meta.id}.bam > ${meta.id}.flagstat.txt
    """
}
