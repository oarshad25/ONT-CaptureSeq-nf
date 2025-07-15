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

    // junction bed file, optional
    path junc_bed

    // preset, -x argument to Minimap2, e.g. "splice"
    val x

    // kmer, -k argument
    val k

    // -u argument
    val u

    // -G argument
    val G

    // -I argument
    val I

    //--cs argument
    val cs

    // any additional options
    val extra_opts

    output:
    tuple val(meta), path("${meta.id}.sam"), emit: sam

    script:
    def preset_flag = x ? "-ax ${x}" : "-ax splice"
    def kmer_flag = k ? "-k ${k}" : ""
    def strand_flag = u ? "-u ${u}" : ""
    def max_intron_length_flag = G ? "-G ${G}" : ""
    def indexing_batchsize_flag = I ? "-I ${I}" : ""
    def cs_flag = cs ? "--cs=${cs}" : ""
    def junc_bed_flag = junc_bed.name != 'NO_FILE' ? "--junc-bed ${junc_bed}" : ''

    """
    minimap2 \\
        -t ${task.cpus} \\
        ${preset_flag} \\
        ${kmer_flag} \\
        ${strand_flag} \\
        ${max_intron_length_flag} \\
        ${indexing_batchsize_flag} \\
        ${cs_flag} \\
        --secondary=yes \\
        ${extra_opts} \\
        ${junc_bed_flag} \\
        ${genome} \\
        ${fastq} > ${meta.id}.sam
    """
}
