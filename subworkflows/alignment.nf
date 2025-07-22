#!/usr/bin/env nextflow

/*
Alignment subworkflow

* Align input reads to input genome.
* Convert alignment ouitput to sam to sorted BAM with samtools and index
  and generate alignment statistics with flagstat
* Optionally calculate read distribution statistics with RSeQC
*/

include { MINIMAP2_INDEX } from "../modules/minimap2_index.nf"
include { MINIMAP2 } from "../modules/minimap2"
include { SAMTOOLS } from "../modules/samtools"

include { RSEQC } from "./rseqc"

workflow ALIGNMENT {
    take:
    reads // [ val(meta), path(fasta) ]
    genome // path to genome fasta
    annotation // path to annotation gtf
    skip_save_minimap2_index // Boolean, skip saving the minimap2 index. This flag controls whether to run indexing seperately
    minimap2_indexing_extra_opts // string, any additional options to pass to indexing process
    minimap2_junc_bed // path to optional junction bed file for MiniMap2
    minimap2_x // params.minimap2_x
    minimap2_k
    minimap2_u
    minimap2_G
    minimap2_I
    minimap2_cs
    minimap2_extra_opts
    skip_rseqc // Boolean, whether to skip RSeQC read distribution calculations

    main:

    /*
    * Index the genome (if flag to prebuild index is set)
    */

    // index the genome
    if (!skip_save_minimap2_index) {
        MINIMAP2_INDEX(
            genome,
            minimap2_x,
            minimap2_k,
            minimap2_I,
            minimap2_indexing_extra_opts,
        )
        ch_minimap2_ref = MINIMAP2_INDEX.out.index
    }
    else {
        ch_minimap2_ref = genome
    }

    // align reads to genome with MiniMap
    MINIMAP2(
        reads,
        ch_minimap2_ref,
        minimap2_junc_bed,
        minimap2_x,
        minimap2_k,
        minimap2_u,
        minimap2_G,
        minimap2_I,
        minimap2_cs,
        minimap2_extra_opts,
    )

    /*
    Use samtools to:
     * convert aligned reads from SAM to BAM
     * sort and index BAM
     * generate alignment statistics
    */

    SAMTOOLS(MINIMAP2.out.sam)

    /*
    * Aligned reads QC
    */

    //initialise empty channel to contain rseqc read distribution output
    rseqc_read_dist_ch = Channel.empty()

    // calculate read distribution of aligned reads with RSeQC if parameter 'skip_resqc' is set
    if (!skip_rseqc) {
        // run RSeQC subworkflow
        RSEQC(SAMTOOLS.out.bambai, annotation)
        // channel with text files containing read distribution calculations
        rseqc_read_dist_ch = RSEQC.out.txt
    }

    emit:
    bambai = SAMTOOLS.out.bambai // sorted and indexed reads: [val(meta), path(bam), path(bai)]
    flagstat = SAMTOOLS.out.flagstat // alignment flagstats [val(meta), path(flagstat_file)]
    rseqc_read_dist = rseqc_read_dist_ch // RSeQC read distribution calculations [val(meta), path(read_dist_file)] or Channel.empty(), if skip_rseqc
}
