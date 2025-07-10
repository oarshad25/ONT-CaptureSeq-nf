// Use Nanocomp to create collated summary of QC metrics for list of
// input files

process NANOCOMP {
    tag "${step}"
    label "low"

    publishDir "${params.outdir}/qc/${step}/nanocomp/", mode: 'copy'

    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/nanocomp:1.24.2--pyhdfd78af_0'
        : 'quay.io/biocontainers/nanocomp:1.24.2--pyhdfd78af_0'}"

    input:
    // list of fastq files
    path inputfiles

    // string. Used in publishDir to set path
    // allows us to put seperate invocations on nanoplot
    // e.g. on raw or filtered data
    // to be published in different directories
    val step

    output:
    path "*.html", emit: html
    path "*.png", emit: png
    path "*NanoStats.txt", emit: stats_txt

    script:
    // use file simpleName for naming in reports
    def names = inputfiles.collect { fn -> "${fn.simpleName}" }.join(' ')

    // set input type option based on whether input files are fastq, BAM or ONT summary file
    if (inputfiles.every { it.name.endsWith('.txt') }) {
        input_type_opt = "--summary"
    }
    else if (inputfiles.every { it.name.contains('.fastq') }) {
        input_type_opt = "--fastq"
    }
    else if (inputfiles.every { it.name.endsWith('.bam') }) {
        input_type_opt = "--bam"
    }
    else {
        error("Please use only *one* of fastq, bam or Nanopore sequencing summary")
    }
    """
    NanoComp \\
        --threads ${task.cpus} \\
        ${input_type_opt} ${inputfiles} \\
        --names ${names}
    """
}
