// Build minimap index from reference genome

// Can use this process to prebuild genome index to speed up the alignemnet process a bit.
// Minimap indexing is fast, building and saving the index leads to even shorter startup time

process MINIMAP2_INDEX {
    label 'medium'

    conda "bioconda::minimap2=2.30 bioconda::samtools=1.23"
    container "${workflow.containerEngine == 'apptainer'
        ? 'oras://community.wave.seqera.io/library/minimap2_samtools:a81ff6397062f3f9'
        : 'community.wave.seqera.io/library/minimap2_samtools:f2e55a1cd407fcaf'}"

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
