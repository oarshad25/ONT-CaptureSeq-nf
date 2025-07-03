// helper function to print help message
def helpMessage() {
    log.info(
        """
        Usage:

        nextflow run main.nf --inputdir <path> [options]

        Required Arguments:

        Input Data:
        --inputdir <path>                   Path to parent read directory containing demultiplexed
                                            barcode subdirectories with fastq reads for each barcode.
        or
        --samplesheet <path>                Path to sample manifest 'csv' file.
                                            Must have columns: id, path.

        Reference Data:
        --genome <path>                     Reference genome FASTA to use for alignment.
        --annotation <path>                 Reference annotation as GTF.

        Optional Arguments:
        --outdir <path>                     Path to output directory.
        --min_reads_per_sample <int>        Threshold to filter samples

        Tool Options:
        --is_fastq_rich <Boolean>           Parameter to set input data type argument for nanoplot
                                            specifying whether input fastq files are in 'rich' format.
        """.stripIndent()
    )
}
