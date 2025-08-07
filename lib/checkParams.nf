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
        error("Please provide either 'read directory' or sample manifest with --inputdir OR --samplesheet")
    }

    // Should specify --inputdir OR --samplesheet, but not both
    if (params.inputdir && params.samplesheet) {
        error("Must specify either --inputdir OR --samplesheet, but not both")
    }

    if (!params.genome) {
        error("No genome specified. Please provide genome fasta with --genome")
    }

    if (!params.annotation) {
        error("No annotation GTF specified. Please provide GTF with --annotation")
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

    // if input directory is provided, check that it is not empty
    if (params.inputdir && file(params.inputdir).listFiles().size() == 0) {
        error("Input directory: `${params.inputdir}` is empty")
    }

    //if restrander is to be run, check config exists

    if (params.run_restrander && !file(params.restrander_config).exists()) {
        error("Specified restrander config: `${params.restrander}` does not exist.")
    }

    // check that no external junction bed file is provided for minimap via 'params.minimap2_junc_bed' if flag to use annotation gtf 'params.alignment_use_annotation' is set
    if (params.alignment_use_annotation && file(params.minimap2_junc_bed).exists() && file(params.minimap2_junc_bed).name != 'NO_FILE') {
        error("Mutually exclusive parameters set. Parameter to use annotation gtf in alignment 'params.alignment_use_annotation' is set to true as well as providing external junction bed file for minimap2 via 'params.minimap2_junc_bed' ${params.minimap2_junc_bed}")
    }

    // check that if minimap2_junc_bonus is provided than either flag 'alignment_use_annotation' is set or junction bed file is provided
    if (params.minimap2_junc_bonus && !params.alignment_use_annotation && file(params.minimap2_junc_bed).name == 'NO_FILE') {
        error("params.minimap2_junc_bonus has been provided. However neither 'params.alignment_use_annotation' is true nor external junction bed file for minimap2 via 'params.minimap2_junc_bed' has been provided")
    }
}
