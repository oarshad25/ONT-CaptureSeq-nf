// helper function to print help message
def helpMessage() {
    log.info(
        """
        Usage:

        nextflow run main.nf -profile <docker/apptainer> --samplesheet <path> --genome <path> --annotation <path>[options]

        Required Arguments:
        -------------------

        Input Data:
        --inputdir <path>                           Path to parent read directory containing demultiplexed
                                                    barcode subdirectories with fastq reads for each barcode.
        or
        --samplesheet <path>                        Path to sample manifest 'csv' file.
                                                    Must have columns: id, path.

        Reference Data:
        --genome <path>                             Reference genome FASTA to use for alignment.
        --annotation <path>                         Reference annotation as GTF.

        Optional Arguments:
        -------------------
        --outdir <path>                             Path to output directory.

        Additional Options:
        -------------------

        QC:
        ------

        --min_reads_per_sample <int>                Threshold to filter samples
        --is_fastq_rich <Boolean>                   Parameter to set input data type argument for nanoplot
                                                    specifying whether input fastq files are in 'rich' format.
        --skip.read_filtering <Boolean>             Skip filtering of reads on length and average quality
        --min_length <integer>                      Minimum read length threshold to filter reads
        --min_qual <float>                          Minimum read average quality threshold to filter reads

        Restrander:
        --run_restrander <Boolean>                  Run restrander
        --restrander_config <path>                  Path to restrander config json

        ALIGNMENT:
        ------

        Minimap2:
        --skip_save_minimap2_index <Boolean>        Whether to skip prebuilding Minimap2 index
        --minimap2_indexing_extra_opts <string>     Any additional options to pass to Minimap2 indexing process
        --minimap2_junc_bed <path>                  Path to optional junction bed file (used by --junc-bed)
        --minimap2_x <string>                       Preset (-x) option e.g. "splice"
        --minimap2_k <integer>                      Kmer size (-k) option e.g. 14
        --minimap2_u <character>                    Strand (-u) option e.g. "f"
        --minimap2_G <string>                       Max intron length (-G) option e.g. "500k"
        --minimap2_I <string>                       Indexing batch size (-I) option e.g. "8G"
        --minimap2_cs <string>                      (-cs) option ("short" or "long") e.g. "long"
        --minimap2_extra_opts <string>              Any additional options to pass to Minimap e.g. "--splice-flank=no"

        RSeQC:
        --skip_rseqc <Boolean>                      Skip read distribution calculation with RSeQC

        Alignment additional:
        --filter_bam_mapped <Boolean>               Whether to filter alignment BAM to mapped reads only

        ISOFORM DISCOVERY:
        ------

        --skip_isoform_discovery <Boolean>          Whether to skip running isoform discovery subworkflow
        --isoform_discovery_method <string>         Method to use for isoform discovery (options: 'isoquant' or 'flair')

        Isoquant:
        --isoquant_complete_genedb <Boolean>        Whether to set option --complete_genedb
        --isoquant_extra_opts <string>              Any additional command line options to pass to IsoQuant

        Flair:
        --flair_collapse_extra_opts <string>        Any additional command line options to pass to flair collapse module
        --flair_align_reads_manifest <string>       Path to sample manifest for flair align module.

        """.stripIndent()
    )
}
