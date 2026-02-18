// Use seqkit stats to generate summary statistics for input fastq
// https://bioinf.shenwei.me/seqkit/usage/#stats
// Also renames input fastq so that first column of output stats file (filename)
// just contains sample id as multiqc picks sample name from this column

process SEQKIT_STATS {
    tag "${meta.id}"
    label "low"

    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    publishDir "${params.outdir}/qc/${step}/seqkit_stats", pattern: '*.tsv', mode: 'copy'

    input:
    // input fastq
    tuple val(meta), path(fastq)

    // string. Used in publishDir to set path
    // allows us to put seperate invocations on nanoplot
    // e.g. on raw or filtered data
    // to be published in different directories
    val step

    output:
    tuple val(meta), path("${meta.id}.${step}_seqkit_stats.tsv"), emit: stats

    script:
    """
    # rename input fastq with just sample id for multiqc (e.g. sample01.full_length.fastq -> sample01.fastq)
    # this ensures file name is just the sample id (as multiqc picks seqkit stats sample name from file name
    # which is the first column in the seqkit stats tsv file)
    mv ${fastq} ${meta.id}.fastq

    seqkit stats \\
        --all \\
        --tabular \\
        --threads ${task.cpus} \\
        --basename \\
        ${meta.id}.fastq > "${meta.id}.${step}_seqkit_stats.tsv"
    """
}
