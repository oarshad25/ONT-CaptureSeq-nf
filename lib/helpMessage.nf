// helper function to print help message
def helpMessage() {
    log.info(
        """
        Usage:

        nextflow run main.nf --inputdir <path> [options]

        Required Arguments:

        Input Data:
        --inputdir <path>           Path to parent read directory containing demultiplexed
                                    barcode subdirectories with fastq reads for each barcode.
        or
        --samplesheet <path>        Path to sample manifest 'csv' file.
                                    Must have columns: id, path

        Reference Data:
        --genome <path>             Reference genome FASTA to use for alignment.
        --annotation <path>         Reference annotation as GTF.

        Optional Arguments:
        --outdir <path>             Path to output directory.
        """.stripIndent()
    )
}
