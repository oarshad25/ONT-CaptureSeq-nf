// unzip input .gz file
// used to decompress reference files if required

process GUNZIP {
    label "low"

    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:24.04'
        : 'quay.io/biocontainers/ubuntu:24.04'}"

    input:
    // path to *.gz to be uncompressed
    path gz_file

    output:
    // uncompressed file without .gz extension
    path "${gz_file.baseName}", emit: gunzip_file

    script:
    """
    gzip -df ${gz_file}
    """
}
