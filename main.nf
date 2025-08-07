#!/usr/bin/env nextflow

// include helper functions
include { helpMessage } from "./lib/helpMessage"
include { checkParams } from "./lib/checkParams"
include { pipelineHeader } from "./lib/pipelineHeader"
include { findFastqFiles } from "./lib/utilities.nf"

// include process modules
include { CONCAT_FASTQ } from "./modules/concat_fastq"
include { FILTER_READS } from "./modules/filter_reads"
include { RESTRANDER } from "./modules/restrander"
include { MULTIQC } from "./modules/multiqc"
include { ISOQUANT } from "./modules/isoquant"
include { ISOQUANT_VISUALIZE } from "./modules/isoquant_visualize"

// include subworkflows
include { PREPARE_REFERENCE_FILES } from "./subworkflows/prepare_reference_files"
include { READ_QC as RAW_READ_QC ; READ_QC as FILTERED_READ_QC ; READ_QC as RESTRANDED_READ_QC } from "./subworkflows/read_qc"
include { ALIGNMENT } from "./subworkflows/alignment"

workflow {

    // Show help message if --help is run
    if (params.help) {
        helpMessage()
        System.exit(1)
    }

    // TODO: validate parameters using nf-schema instead of using own function
    // https://nextflow-io.github.io/nf-schema/latest/parameters/validation/
    checkParams()

    // print pipeline header
    pipelineHeader()

    /*
    * CREATE INPUT CHANNELS
    */

    // If the --inputdir option was provided
    if (params.inputdir) {

        log.info("USING INPUT DIRECTORY: ${params.inputdir}")

        /*
        make input data channel with the sample metamap and path to fastq
        directory of reads
        [meta, fastqdir]
        Take 'params.inputdir' and create a tuple channel [meta, fastqdir]
        - meta: a metadata map. meta.id is he sample id which is taken from
          the name of each sample subdirectory in params.inputdir e.g. [id: barcode01]
        - fastqdir: path to each sample subdirectory containg the reads for the sample
        From inputdir, we get the list of subdirectories, which is flattened
        to emit an individual entry for each item. Each subdirectory is then
        mapped to a sample id which is simply the basename of the subdirectory
        which is stored in the meta map which is the first element of the channel.
        Path to the sample subdirectory containing the reads for the sample
        is the second element of the channel
        Unclassified reads subdirectory is filtered out
        */
        Channel.fromPath(params.inputdir, checkIfExists: true)
            .map { inputdir -> file(inputdir).listFiles() }
            .flatten()
            .map { fastqdir ->
                def id = fastqdir.Name
                def meta = [id: id]
                [meta, fastqdir]
            }
            .filter { meta, _fastqdir -> !meta.id.contains("unclassified") }
            .set { input_ch }
    }
    else {
        // otherwise samplesheet was provided
        log.info("USING SAMPLESHEET: ${params.samplesheet}")

        // create input channel from sample sheet

        // TODO: Validate samplesheet
        // Add validation of samplesheet to ensure it has requisite columns

        // Make a channel which includes
        // The sample name from the 'id' column
        // and path to fastq directory for reads from 'fastqdir' column
        // [meta, fastqdir]
        Channel.fromPath(params.samplesheet, checkIfExists: true)
            .splitCsv(header: true)
            .flatten()
            .map { row -> [row.id, file(row.fastqdir)] }
            .map { id, fastqdir ->
                def meta = [id: id]
                [meta, fastqdir]
            }
            .set { input_ch }
    }

    /*
    * Setup reference channels
    */

    // replace reference channel creation below with subworkflow PREPARE_REFERENCE_FILES
    // create reference genome fasta channel
    //genome_ch = file(params.genome, checkIfExists: true)

    // create reference genome fasta channel
    //annotation_ch = file(params.annotation, checkIfExists: true)

    PREPARE_REFERENCE_FILES(params.genome, params.annotation)

    // reference genome fasta channel
    genome_ch = PREPARE_REFERENCE_FILES.out.genome_ch
    // reference annotation channel
    annotation_ch = PREPARE_REFERENCE_FILES.out.annotation_ch

    /*
    * CONCATENATE FQS PER SAMPLE
    */

    // create a channel colecting all fastq files in the sample directory into a list
    // [meta, fastq_files]

    // The filter operator:
    // Filter out any samples in input_ch for which no fastq files found,
    // preventing channel with empty components being passed on

    input_ch
        .map { meta, fastqdir ->
            // get fastq files in fastqdir
            def fastq_files = findFastqFiles(fastqdir)
            //output a tuple of metadata map and list of fastq_files for sample
            tuple(meta, fastq_files)
        }
        .filter { it -> it[1] }
        .set { fastq_files_ch }

    // concatenate fastq files per sample
    CONCAT_FASTQ(fastq_files_ch)

    /*
    * REMOVE SAMPLES WITH FEW READS
    */

    // filter out samples in the channel with less than `params.min_reads_per_sample` reads
    // this is mainly to get rid of unassigned barcodes

    // channel containing the merged fastq for samples that pass read count filter
    // [meta, fastq]
    reads_ch = CONCAT_FASTQ.out.merged_reads
        .map { meta, fastq -> [meta, fastq, fastq.countFastq()] }
        .filter { _meta, _fastq, numreads -> numreads > params.min_reads_per_sample }
        .map { meta, fastq, _numreads -> tuple(meta, fastq) }


    /*
    * QC on raw reads
    */
    RAW_READ_QC(reads_ch, "raw", params.is_fastq_rich)

    /*
    * Filter raw reads
    */

    // initialise to upstream channel
    filtered_reads_ch = reads_ch

    // filter raw reads on length and quality
    if (!params.skip_read_filtering) {
        FILTER_READS(reads_ch, params.min_length, params.min_qual)

        filtered_reads_ch = FILTER_READS.out.filtered_reads

        // QC stats on filtered reads
        FILTERED_READ_QC(filtered_reads_ch, "filtered", params.is_fastq_rich)
    }

    /*
    * Restrand filtered reads
    */

    // initialise restranded reads
    full_length_reads_ch = filtered_reads_ch

    if (params.run_restrander) {

        // channel containing path to restrander config file
        restrander_config_ch = file(params.restrander_config, checkIfExists: true)

        // reorient reads with restrander
        RESTRANDER(filtered_reads_ch, restrander_config_ch)

        full_length_reads_ch = RESTRANDER.out.full_length_reads

        // QC on restranded reads
        RESTRANDED_READ_QC(full_length_reads_ch, "full_length", params.is_fastq_rich)
    }

    // remove any samples that after preprocessing have less than 'min_reads_per_sample'
    full_length_reads_ch
        .map { meta, fastq -> [meta, fastq, fastq.countFastq()] }
        .filter { _meta, _fastq, numreads -> numreads > params.min_reads_per_sample }
        .map { meta, fastq, _numreads -> tuple(meta, fastq) }
        .set { processed_reads_ch }

    /*
    * ALIGNMENT
    */

    minimap2_junc_bed_ch = file(params.minimap2_junc_bed, checkIfExists: true)

    // align reads to genome with MiniMap and generate read QC statistics
    ALIGNMENT(
        processed_reads_ch,
        genome_ch,
        annotation_ch,
        params.alignment_use_annotation,
        params.skip_save_minimap2_index,
        params.minimap2_indexing_extra_opts,
        minimap2_junc_bed_ch,
        params.minimap2_x,
        params.minimap2_k,
        params.minimap2_u,
        params.minimap2_G,
        params.minimap2_I,
        params.minimap2_cs,
        params.minimap2_junc_bonus,
        params.minimap2_extra_opts,
        params.skip_rseqc,
        params.filter_bam_mapped,
    )

    // aligned reads, sorted and indexed bam: [val(meta), path(bam), path(bai)]
    aligned_reads_ch = ALIGNMENT.out.bambai

    /*
    * Alignment MultiQC report
    */

    // alignment flagstats input channel for MultiQC
    multiqc_flagstat_ch = ALIGNMENT.out.flagstat.collect { it -> it[1] }

    // nanostats input channel for multiQC
    multiqc_alignment_nanostats_ch = ALIGNMENT.out.nanostats.collect { it -> it[1] }

    // channel with text files containing read distribution calculations
    multiqc_rseqc_ch = ALIGNMENT.out.rseqc_read_dist.collect { it -> it[1] }.ifEmpty([])

    // combine samtools flagstat, nanostats and RSeQC read distribution files
    multiqc_alignment_input_files_ch = multiqc_flagstat_ch
        .mix(
            multiqc_alignment_nanostats_ch,
            multiqc_rseqc_ch,
        )
        .collect()

    // run MultiQC on aligned read statistics
    MULTIQC(multiqc_alignment_input_files_ch, "aligned")

    /*
    * ISOQUANT
    */

    ISOQUANT(
        aligned_reads_ch.collect { it -> it[1] },
        aligned_reads_ch.collect { it -> it[2] },
        genome_ch,
        annotation_ch,
        params.isoquant_complete_genedb,
        params.isoquant_extra_opts,
    )

    // Note: Following section is a stub as IsoQuant visualisation tool does not work
    // TODO: look at GitHub issues to see if IsoQuant visualisation tool fixed
    if (params.isoquant_visualize_genelist) {
        isoquant_visualize_genelist_ch = Channel.fromPath(params.isoquant_visualize_genelist, checkIfExists: true)
        ISOQUANT_VISUALIZE(
            ISOQUANT.out.output_dir,
            isoquant_visualize_genelist_ch,
            annotation_ch,
        )
    }

    /*
    * Workflow event handlers
    */

    if (params.workflow_event_handlers) {
        // Print workflow execution summary
        workflow.onComplete = {
            println("")
            println("=======================================================================================")
            println("Workflow execution summary")
            println("=======================================================================================")

            if (workflow.success) {
                println("Command line   : ${workflow.commandLine}")
                println("Completed at   : ${workflow.complete}")
                println("Duration       : ${workflow.duration}")
                println("Work directory : ${workflow.workDir}")
                println("\n")
                println("Parameters")
                println("==========")
                params.each { k, v ->
                    println("${k}: ${v}")
                }
            }
            else {
                println("Exit status    : ${workflow.exitStatus}")
            }
        }

        workflow.onError = {
            println("Pipeline execution stopped with the following message: ${workflow.errorMessage}")
        }
    }
}
