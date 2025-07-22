// Build minimap index from reference genome

// Can use this process to prebuild genome index to speed up the alignemnet process a bit.
// Minimap indexing is fast, building and saving the index leads to even shorter startup time

process MINIMAP2_INDEX {
    label 'medium'

    conda "bioconda::minimap2=2.28"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/minimap2:2.28--h577a1d6_4'
        : 'quay.io/biocontainers/minimap2:2.28--h577a1d6_4'}"

    input:
    // genome fasta file
    path genome

    // preset, -x argument to Minimap2, e.g. "splice"
    val x

    // kmer, -k argument
    val k

    // -I argument
    val I

    // any additional options
    val extra_opts

    output:
    path "*.mmi", emit: index

    script:
    def preset_flag = x ? "-x ${x}" : "-x splice"
    def kmer_flag = k ? "-k ${k}" : ""
    def indexing_batchsize_flag = I ? "-I ${I}" : ""

    """
    minimap2 \\
        -t ${task.cpus} \\
        -d ${genome.baseName}.mmi \\
        ${preset_flag} \\
        ${kmer_flag} \\
        ${indexing_batchsize_flag} \\
        ${extra_opts} \\
        ${genome}
    """
}
