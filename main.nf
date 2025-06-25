#!/usr/bin/env nextflow

// include helper functions
include { helpMessage } from "./lib/helpMessage"
include { checkParams } from "./lib/checkParams"
include { pipelineHeader } from "./lib/pipelineHeader"

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
    ============================================================================
    Create input channels
    ============================================================================
    */

    // If the --inputdir option was provided
    if (params.inputdir) {

        log.info("USING INPUT DIRECTORY: ${params.inputdir}")

        /*
        set input data channel [meta, fastqdir]
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

        input_ch.view()
    }
    else {
        // otherwise samplesheet was provided
        log.info("USING SAMPLESHEET: ${params.samplesheet}")

        // create input channel from sample sheet [meta, fastqdir]

        // TODO: Validate samplesheet
        // Add validation of samplesheet to ensure it has requisite columns

        // Make a channel which includes
        // The sample name from the 'id' column
        // and path to fastq directory for reads from 'fastqdir' column
        Channel.fromPath(params.samplesheet, checkIfExists: true)
            .splitCsv(header: true)
            .flatten()
            .map { row -> [row.id, file(row.fastqdir)] }
            .map { id, fastqdir ->
                def meta = [id: id]
                [meta, fastqdir]
            }
            .set { input_ch }

        input_ch.view()
    }

    // create reference genome fasta channel
    genome_ch = file(params.genome, checkIfExists: true)

    // create reference genome fasta channel
    annotation_ch = file(params.annotation, checkIfExists: true)
}
