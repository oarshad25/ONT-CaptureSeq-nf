// helper function to print help message
def helpMessage() {
    log.info(
        """
        Usage:

        nextflow run main.nf --inputdir <path> [options]

        Required Arguments:

        Input Data:
        --inputdir              Path to parent read directory containing demultiplexed
                                barcode subdirectories with fastq reads for each barcode.
        or
        --samplesheet           Path to sample manifest 'csv' file.
                                Must have columns: id, path

        Reference Data:
        --genome                Reference genome FASTA to use for alignment.
        --annotation            Reference annotation as GTF.

        Optional Arguments:
        --outdir                Path to output directory.
        """.stripIndent()
    )
}
