#!/usr/bin/env nextflow

// Generate read QC statistics

include { NANOPLOT } from "../modules/nanoplot"
include { NANOCOMP } from "../modules/nanocomp"
include { MULTIQC } from "../modules/multiqc"

workflow READ_QC {
    take:
    reads_ch // tuple val(meta), path(read); queue channel with tuples of metamap and fastq/bam
    step // string. Allows us to publish seperate invocations of workflow for e.g on raw or filtered data, in different directories
    is_fastq_rich // Boolean. Used in nanoplot process to set input file type argument, whether input fastq is in 'rich' format

    main:
    NANOPLOT(reads_ch, step, is_fastq_rich)

    // sort the reads channel by sample id
    // sorting is done so that nanocomp report is sorted by sample id
    reads_ch
        .toSortedList { tup1, tup2 -> tup1[0].id <=> tup2[0].id }
        .flatMap()
        .set { sorted_reads_ch }

    NANOCOMP(sorted_reads_ch.collect { it -> it[1] }, step)

    def multiqc_input_ch = NANOPLOT.out.txt.collect { it -> it[1] }

    // dummy multiqc config file channel
    Channel.fromPath("${projectDir}/assets/NO_FILE").set { multiqc_config_ch }

    MULTIQC(multiqc_input_ch, multiqc_config_ch, step)

    emit:
    nanostats = NANOPLOT.out.txt // queue channel of sample nanostat txt files: [val(meta), path(nanostats_files)] file(sample1_NanoStats.txt)
}
