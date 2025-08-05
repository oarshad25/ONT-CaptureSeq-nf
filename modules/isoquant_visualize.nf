/*
Use IsoQuant's visualization tool to explore transcript usage and splicing
patterns for genes of interest
*/

process ISOQUANT_VISUALIZE {
    label 'low'

    publishDir "${params.outdir}", mode: 'link'

    conda "bioconda::isoquant=3.7.0"
    container "quay.io/biocontainers/isoquant:3.7.0--hdfd78af_0"

    input:
    // directory containing IsoQuant output files
    path isoquant_output_dir

    // Path to a .txt file containing a list of genes, each on its own line.
    path gene_list

    // annotation gtf
    path annotation

    output:
    // directory to save visualisation output files
    path "isoquant_visualize"

    script:
    """
    python /usr/local/bin/visualize.py ${isoquant_output_dir} \\
        --gene_list ${gene_list} \\
        --gtf ${annotation}
    """
}
