// Build minimap index from reference genome

process MINIMAP2 {
    tag "${meta.id}"
    label 'medium'

    conda "bioconda::minimap2=2.28"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/minimap2:2.28--h577a1d6_4'
        : 'quay.io/biocontainers/minimap2:2.28--h577a1d6_4'}"

    input:
    // input reads to align
    tuple val(meta), path(fastq)

    // genome fasta file
    path genome

    // minimap2 extra command line options
    val minimap2_opts

    // junction bed file, optional
    path junc_bed

    output:
    tuple val(meta), path("${meta.id}.sam"), emit: sam

    script:
    def junc_bed_line = junc_bed.name != 'NO_FILE' ? "--junc-bed ${junc_bed}" : ''
    """
    minimap2 \\
        -t ${task.cpus} \\
        ${minimap2_opts} \\
        ${junc_bed_line} \\
        ${genome} \\
        ${fastq} > ${meta.id}.sam
    """
}
