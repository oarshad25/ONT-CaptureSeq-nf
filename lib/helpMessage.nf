// helper function to print help message
def helpMessage() {
    log.info(
        """
        Usage:

        nextflow run main.nf --inputdir <path> [options]

        Required Arguments:
        -------------------

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
        -------------------
        --outdir <path>                     Path to output directory.

        Additional Options:
        -------------------

        QC:
        --min_reads_per_sample <int>        Threshold to filter samples
        --is_fastq_rich <Boolean>           Parameter to set input data type argument for nanoplot
                                            specifying whether input fastq files are in 'rich' format.
        --skip.read_filtering <Boolean>     Skip filtering of reads on length and average quality
        --min_length <integer>              Minimum read length threshold to filter reads
        --min_qual <float>                  Minimum read average quality threshold to filter reads

        Restrander:
        --run_restrander <Boolean>          Run restrander
        --restrander_config <path>          Path to restrander config json
        """.stripIndent()
    )
}
