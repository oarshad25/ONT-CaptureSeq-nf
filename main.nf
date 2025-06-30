#!/usr/bin/env nextflow

// include helper functions
include { helpMessage } from "./lib/helpMessage"
include { checkParams } from "./lib/checkParams"
include { pipelineHeader } from "./lib/pipelineHeader"

// include process modules
include { CONCAT_FASTQ } from "./modules/concat_fastq"

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

        input_ch.view()
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

    // create reference genome fasta channel
    genome_ch = file(params.genome, checkIfExists: true)

    // create reference genome fasta channel
    annotation_ch = file(params.annotation, checkIfExists: true)

    /*
    * CONCATENATE FQS PER SAMPLE
    */

    // https://github.com/asadprodhan/How-to-channel-sequencing-reads-from-multiple-subdirectories-into-nextflow-pipeline
    // create a channel colecting all fastq files in the sample directory into a list
    // [meta, fastq_files]

    fastq_extns = ['.fq', '.fq.gz', '.fastq', '.fastq.gz']

    input_ch
        .map { meta, fastqdir ->
            // list all files in fastqdir
            def all_files = file(fastqdir).listFiles()
            // get files with extensions in fastq_extns
            def fastq_files = all_files.findAll { fn -> fastq_extns.find { fn.name.endsWith(it) } }
            //output a tuple of metadata map and list of fastq_files for sample
            tuple(meta, fastq_files)
        }
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
}
