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

    // If the --inputdir option was provided
    if (params.inputdir) {
    }
    else {
    }
}
