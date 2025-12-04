/*
* use UCSC's genePredToBed to convert genepred to bed
*/

process UCSC_GENEPRED_TO_BED {
    label "single"

    conda "bioconda ucsc-genepredtobed=469"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/ucsc-genepredtobed:469--h664eb37_1'
        : 'quay.io/biocontainers/ucsc-genepredtobed:469--h664eb37_1'}"

    input:
    path genepred

    output:
    path "*.bed", emit: bed

    script:
    """
    genePredToBed \\
        ${genepred}  \\
        ${genepred.simpleName}.bed
    """
}
