// generate QC reports with MultiQC

process MULTIQC {
    tag "${step}"
    label "low"

    publishDir "${params.outdir}/qc/${step}/multiqc/", mode: 'copy'

    container "quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0"

    input:
    path '*'

    // string. Used in publishDir to set path
    // allows us to put seperate invocations on MultiQC
    // e.g. on raw or filtered data
    // to be published in different directories
    val step

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc \\
        --force \\
        .
    """
}
