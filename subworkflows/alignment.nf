#!/usr/bin/env nextflow

/*
Alignment subworkflow

* Align input reads to input genome optionally using input annotation or junction bed file
* Generate sorted/indexed BAM and flagstat with samtools
* QC of aligned reads with cramino
* Optionally filter alignments to primary only (remove secondary, supplementary and unmapped reads)
* Optionally calculate read statistics with RSeQC
*/

include { MINIMAP2_INDEX } from "../modules/minimap2_index"
include { MINIMAP2_PAFTOOLS_GFF2BED } from "../modules/minimap2_paftools_gff2bed"
include { MINIMAP2_SAMTOOLS } from "../modules/minimap2_samtools"
include { CRAMINO } from "../modules/cramino"
include { FILTER_ALIGNMENTS } from "../modules/filter_alignments"

include { RSEQC } from "./rseqc"

workflow ALIGNMENT {
    take:
    reads // [ val(meta), path(fasta) ]
    genome // path to genome fasta
    annotation // path to annotation gtf used by Minimap2
    rseqc_bed // path to annotation in bed format to be used by RSeQC
    rseqc_housekeeping_bed // path to housekeeping gene bed for RSeQC gene body coverage module
    alignment_use_annotation // Boolean, whether to use gene annotation as input to Minimap2 to prioritise on annotated splice junctions
    skip_save_minimap2_index // Boolean, skip saving the minimap2 index. This flag controls whether to run indexing seperately
    minimap2_indexing_extra_opts // string, any additional options to pass to indexing process
    minimap2_junc_bed // path to optional junction bed file for MiniMap2
    minimap2_x // params.minimap2_x
    minimap2_k
    minimap2_u
    minimap2_G
    minimap2_I
    minimap2_cs
    minimap2_junc_bonus
    minimap2_extra_opts
    skip_rseqc // Boolean, whether to skip RSeQC read distribution calculations
    filter_alignments // Boolean, whether to filter alignment BAMs to remove secondary, supplementary and unmapped reads

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

    /*
    * Align reads to the genome with MiniMap2 and generate sorted/indexed BAM + flagstats
    */

    // if flag to use input annotation in alignment is set, convert annotation to bed using Minimap's utility
    // and use this bed file in Minimap alignment to prioritise on annotated spliced junctions
    if (alignment_use_annotation) {
        MINIMAP2_PAFTOOLS_GFF2BED(annotation)
        minimap2_junc_bed = MINIMAP2_PAFTOOLS_GFF2BED.out.bed12
    }

    MINIMAP2_SAMTOOLS(
        reads,
        ch_minimap2_ref,
        minimap2_junc_bed,
        minimap2_x,
        minimap2_k,
        minimap2_u,
        minimap2_G,
        minimap2_I,
        minimap2_cs,
        minimap2_junc_bonus,
        minimap2_extra_opts,
    )

    /*
    * Use cramino for quality assessment of unfiltered BAM files
    */

    CRAMINO(MINIMAP2_SAMTOOLS.out.bambai)


    /*
    * Optionally filter out secondary, supplementray and unmapped reads
    */

    if (filter_alignments) {
        FILTER_ALIGNMENTS(MINIMAP2_SAMTOOLS.out.bambai)
        bambai_ch = FILTER_ALIGNMENTS.out.filtered_bambai
    }
    else {
        bambai_ch = MINIMAP2_SAMTOOLS.out.bambai
    }

    /*
    * Use RSeQC for alignment QC if flag skip_rseqc is not set
    */

    //initialise empty channel to contain rseqc module outputs
    rseqc_read_gc_xls_ch = Channel.empty()
    rseqc_read_dist_ch = Channel.empty()
    rseqc_infer_exp_ch = Channel.empty()
    rseqc_genebody_coverage_txt_ch = Channel.empty()

    // calculate read distribution of filtered aligned reads with RSeQC if alignment filtering is on and parameter 'skip_resqc' is set
    if (filter_alignments && !skip_rseqc) {
        // run RSeQC subworkflow
        RSEQC(bambai_ch, rseqc_bed, rseqc_housekeeping_bed)
        // channel with xls files containing read GC content calculations
        rseqc_read_gc_xls_ch = RSEQC.out.read_gc_xls
        // channel with text files containing read distribution calculations
        rseqc_read_dist_ch = RSEQC.out.read_distribution_txt
        // channel with txt file containing results of infer_experiment module
        rseqc_infer_exp_ch = RSEQC.out.infer_experiment_txt
        // channel with txt file containing results of gene body coverage module
        rseqc_genebody_coverage_txt_ch = RSEQC.out.rseqc_genebody_coverage_txt
    }

    emit:
    bambai = bambai_ch // sorted and indexed reads (optionally filtered to mapped only): [val(meta), path(bam), path(bai)]
    flagstat = MINIMAP2_SAMTOOLS.out.flagstat // alignment flagstats [val(meta), path(flagstat_file)]
    cramino_stats = CRAMINO.out.stats_txt // cramino bam stats [val(meta), path(cramino_stat_file)]
    rseqc_read_gc_xls = rseqc_read_gc_xls_ch // RSeQC read GC content calculation xls file [val(meta), path(read_gc_xls)] or Channel.empty(), if skip_rseqc
    rseqc_read_dist = rseqc_read_dist_ch // RSeQC read distribution calculations [val(meta), path(read_dist_file)] or Channel.empty(), if skip_rseqc
    rseqc_infer_exp = rseqc_infer_exp_ch // RSeQC infer experiment module result [val(meta) path(infer_exp_txt)] or Channel.empty(), if skip_rseqc
    rseqc_genebody_coverage_txt = rseqc_genebody_coverage_txt_ch // RSeQC gene body coverage module result [val(meta) path(genebody_coverage_txt)] or Channel.empty(), if skip_rseqc or no housekeeping gene bed provided
}
