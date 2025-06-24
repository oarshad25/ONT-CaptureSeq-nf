// helper function to print pipeline header
def pipelineHeader() {
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
}
