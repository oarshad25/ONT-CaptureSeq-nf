// helper function to check if required parameters have been provided

include { helpMessage } from "./helpMessage"

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
