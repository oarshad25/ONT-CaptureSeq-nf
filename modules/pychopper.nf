process PYCHOPPER {
    debug true

    tag "${meta.id}"
    label "medium"

    publishDir "${params.outdir}/qc/pychopper/stats", pattern: '*.tsv', mode: 'copy'
    publishDir "${params.outdir}/qc/pychopper/report", pattern: '*_report.pdf', mode: 'copy'
    publishDir "${params.outdir}/fastq/pychopper", pattern: '*.fastq', mode: 'link'

    input:
    // input reads to preprocess
    tuple val(meta), path(fastq)

    // kit name, -k option to pychopper (optional, if primer fasta provided)
    val cdna_kit

    // primer fasta, -b option to pychopper (optional, if cdna_kit provided)
    path primer_fasta

    // pychopper backend, -m option to pychopper
    val backend

    // any additional command line options
    val extra_opts

    output:
    tuple val(meta), path("${meta.id}.full_length.fastq"), emit: full_length_reads
    tuple val(meta), path("${meta.id}.unclassified.fastq"), emit: unclassified_reads
    tuple val(meta), path("${meta.id}.rescued.fastq"), emit: rescued_reads
    tuple val(meta), path("${meta.id}_report.pdf"), emit: pdf
    tuple val(meta), path("${meta.id}.tsv"), emit: stats

    script:
    def cdna_kit_flag = cdna_kit ? "-k ${cdna_kit}" : ""
    def primer_fasta_flag = primer_fasta.name != 'NO_FILE' ? "-b ${primer_fasta}" : ""

    """
    pychopper \\
        -t ${task.cpus} \\
        ${cdna_kit_flag} \\
        ${primer_fasta_flag} \\
        -m ${backend} \\
        -S ${meta.id}.tsv \\
        -r ${meta.id}_report.pdf \\
        -u ${meta.id}.unclassified.fastq \\
        -w ${meta.id}.rescued.fastq \\
        ${extra_opts} \\
        ${fastq} \\
        ${meta.id}.full_length.fastq
    """
}
