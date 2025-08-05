// Generate quality metrics for the input data with nanoplot

process NANOPLOT {
    tag "${meta.id}"
    label "low"

    publishDir "${params.outdir}/qc/${step}/nanoplot/${meta.id}/", mode: 'copy'

    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/nanoplot:1.43.0--pyhdfd78af_1'
        : 'quay.io/biocontainers/nanoplot:1.43.0--pyhdfd78af_0'}"

    input:
    // input file is a .fastq(gz) file, *.bam or *.txt ONT summary file
    tuple val(meta), path(inputfile)

    // string. Used in publishDir to set path
    // allows us to put seperate invocations on nanoplot
    // e.g. on raw or filtered data
    // to be published in different directories
    val step

    // Boolean. Is input fastq in 'rich' format?
    val is_fastq_rich

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.log"), emit: log

    script:
    /*
    * Set input data type option

    * If input file is a text file, assume summary file provided.
    * If it is a BAM file, set input data source flag to be --bam
    * If it is a fastq file ending in (*.fastq(.gz)),
    * set input data source based on Boolean parameter 'nanoplot_is_fastq_rich'
    */
    // if summary file provided
    if (inputfile.name.endsWith(".txt")) {
        input_type_opt = "--summary"
    }
    else if (inputfile.name.endsWith(".fastq") || inputfile.name.endsWith(".fastq.gz")) {
        // check if is_fastq_rich is set
        if (!is_fastq_rich) {
            input_type_opt = "--fastq"
        }
        else {
            input_type_opt = "--fastq_rich"
        }
    }
    else if (inputfile.name.endsWith(".bam")) {
        input_type_opt = "--bam"
    }
    else {
        error("Unrecognised input file: ${inputfile}")
    }
    """
    NanoPlot \\
        --threads ${task.cpus} \\
        --no_static \\
        --prefix "${meta.id}_" \\
        ${input_type_opt} ${inputfile}
    """
}
