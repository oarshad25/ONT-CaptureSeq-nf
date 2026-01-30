// filter reads based on length and quality using seqkit seq

process FILTER_READS {
    tag "${meta.id}"
    label "low"

    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    publishDir "${params.outdir}/fastq/filtered", pattern: '*.filtered.fastq', mode: 'link'

    input:
    // input reads to filter
    tuple val(meta), path(fastq)

    // filter out reads shorter than min_lenth
    val min_length

    // filter out reads longer than max_length
    val max_length

    // filter out reads with an average quality below min_qual
    val min_qual

    output:
    tuple val(meta), path("*.filtered.fastq"), emit: filtered_reads

    script:
    """
    seqkit seq ${fastq} \\
        --out-file "${meta.id}.filtered.fastq" \\
        --min-len ${min_length} \\
        --max-len ${max_length} \\
        --min-qual ${min_qual}
    """
}
