#!/usr/bin/env nextflow

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

// helper function to check if required parameters have been provided
def checkParams() {
    // must specify one of --inputdir or --samplesheet
    if (!params.inputdir && !params.samplesheet) {
        log.info(
            """
            "Please provide either 'read directory' or sample manifest with --inputdir OR --samplesheet"
            """.stripIndent()
        )
        helpMessage()
        System.exit(1)
    }

    // Should specify --inputdir OR --samplesheet, but not both
    if (params.inputdir && params.samplesheet) {
        log.info(
            """
            Must specify either --inputdir OR --samplesheet, but not both
            """.stripIndent()
        )
        helpMessage()
        System.exit(1)
    }

    if (!params.genome) {
        log.info(
            """
            "No genome specified. Please provide genome fasta with --genome"
            """.stripIndent()
        )
        helpMessage()
        System.exit(1)
    }

    if (!params.annotation) {
        log.info(
            """
            "No annotation GTF specified. Please provide GTF with --annotation"
            """.stripIndent()
        )
        helpMessage()
        System.exit(1)
    }
}

workflow {

    // Show help message if --help is run
    if (params.help) {
        helpMessage()
        System.exit(1)
    }

    // TODO: validate parameters using nf-schema instead of oneself
    // https://nextflow-io.github.io/nf-schema/latest/parameters/validation/
    checkParams()

    // Print pipeline header
    // TODO: Use text to ascii art generator for pipeline name in header
    // Set a header made using https://patorjk.com/software/taag
    // e.g. https://github.com/vibbits/nextflow-workshop/blob/main/project-solutions/main.nf
    log.info(
        """
        =======================================================================================
        ONT-CaptureSeq-nf
        =======================================================================================

        =======================================================================================
        Workflow run parameters
        =======================================================================================
        Read Directory              : ${params.inputdir}
        Sample Manifest             : ${params.samplesheet}
        Genome FASTA                : ${params.genome}
        Genome Annotation GTF       : ${params.annotation}
        Results Directory           : ${params.outdir}
        =======================================================================================
        """.stripIndent()
    )

    // If the --inputdir option was provided
    if (params.inputdir) {
    }
    else {
    }
}
