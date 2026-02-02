// get full length reads with restrander

process RESTRANDER {
    tag "${meta.id}"
    label "low"

    container 'genomicpariscentre/restrander:1.1.1'

    publishDir "${params.outdir}/qc/restranded/restrander_stats", pattern: '*.json', mode: 'copy'
    publishDir "${params.outdir}/fastq/restranded/", pattern: '*.restranded.fastq', mode: 'link'

    input:
    // input reads to restrand
    tuple val(meta), path(fastq)

    // config file
    path config

    output:
    tuple val(meta), path("${meta.id}.restranded.fastq"), emit: restranded_reads
    tuple val(meta), path("${meta.id}.stats.json"), emit: restrander_stats

    script:
    """
    restrander \\
        ${fastq} \\
        ${meta.id}.restranded.fastq \\
        ${config} \\
        > ${meta.id}.stats.json
    """
}
