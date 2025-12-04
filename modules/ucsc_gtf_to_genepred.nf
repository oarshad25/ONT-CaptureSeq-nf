/*
* use UCSC's gtfToGenePred to convert a GTF file to genePred
*/

process UCSC_GTF_TO_GENEPRED {
    label "single"

    conda "bioconda ucsc-gtftogenepred=469"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/ucsc-gtftogenepred:469--h664eb37_1'
        : 'quay.io/biocontainers/ucsc-gtftogenepred:469--h664eb37_1'}"

    input:
    path gtf

    output:
    path "*.genepred", emit: genepred

    script:
    """
    gtfToGenePred \\
        ${gtf}  \\
        ${gtf.simpleName}.genepred
    """
}
