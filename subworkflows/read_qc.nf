#!/usr/bin/env nextflow

// Generate read QC statistics

include { NANOPLOT } from "../modules/nanoplot"
include { NANOCOMP } from "../modules/nanocomp"
include { MULTIQC } from "../modules/multiqc"

workflow READ_QC {
    take:
    reads_ch // tuple val(meta), path(fastq); queue channel with tuples of metamap and fastq
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

    MULTIQC(multiqc_input_ch, step)

    emit:
    nanostats_files = NANOPLOT.out.txt.map { it -> it[1] } // queue channel of sample nanostat txt files file(sample1_NanoStats.txt)
}
