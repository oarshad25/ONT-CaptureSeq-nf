// helper function to check if required parameters have been provided

include { helpMessage } from "./helpMessage"

def checkParams() {

    /*
    validate parameters
    -------------------
    check if required parameters have been provided
    */

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

    /*
    validate paths
    -------------------
    check if provided input paths for input directory or samplesheet and
    references exist
    */

    //check if provided input paths exist

    if (params.inputdir && !file(params.inputdir).exists()) {
        error("Input directory: `${params.inputdir}` does not exist")
    }

    if (params.samplesheet && !file(params.samplesheet).exists()) {
        error("Sample manifest: `${params.samplesheet}` does not exist")
    }

    if (params.genome && !file(params.genome).exists()) {
        error("Genome: `${params.genome}` does not exist")
    }

    if (params.annotation && !file(params.annotation).exists()) {
        error("Annotation: `${params.annotation}` does not exist")
    }
}
