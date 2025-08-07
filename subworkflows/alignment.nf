#!/usr/bin/env nextflow

/*
Alignment subworkflow

* Align input reads to input genome optionally using input annotation or junction bed file
* Convert alignment output to sam to sorted BAM with samtools and index
  and generate alignment statistics with flagstat
* Generate aligned reads length and quality histograms with seqkit bam
* Optionally calculate read distribution statistics with RSeQC
* Optionally filter alignments to mapped reads only
*/

include { MINIMAP2_INDEX } from "../modules/minimap2_index"
include { MINIMAP2_PAFTOOLS_GFF2BED } from "../modules/minimap2_paftools_gff2bed"
include { MINIMAP2 } from "../modules/minimap2"
include { SAMTOOLS } from "../modules/samtools"
include { CRAMINO } from "../modules/cramino"
include { NANOPLOT } from "../modules/nanoplot"
include { FILTER_BAM_MAPPED } from "../modules/filter_bam_mapped"

include { RSEQC } from "./rseqc"

workflow ALIGNMENT {
    take:
    reads // [ val(meta), path(fasta) ]
    genome // path to genome fasta
    annotation // path to annotation gtf used by RSeqC and optionally Minimap2
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
    filter_bam_mapped // Boolean, whether to get mapped reads only from aligned BAM

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
    * Align reads to the genome with MiniMap2
    */

    // if flag to use input annotation in alignment is set, convert annotation to bed using Minimap's utility
    // and use this bed file in Minimap alignment to prioritise on annotated spliced junctions
    if (alignment_use_annotation) {
        MINIMAP2_PAFTOOLS_GFF2BED(annotation)
        minimap2_junc_bed = MINIMAP2_PAFTOOLS_GFF2BED.out.bed12
    }

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
        minimap2_junc_bonus,
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
    * Use cramino for quality assessment of BAM files
    */

    // Alignment QC statistics are computed on unfiltered BAMs

    CRAMINO(SAMTOOLS.out.bambai)

    /*
    * Run nanoplot to get metrics on alignment bams
    */

    def nanoplot_input_ch = SAMTOOLS.out.bambai.map { meta, bam, _bai -> tuple(meta, bam) }
    NANOPLOT(nanoplot_input_ch, "aligned", "")

    /*
    * Use RSeQC to generate read distribution
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

    /*
    * Optionally filter out unmapped reads
    */

    if (filter_bam_mapped) {
        FILTER_BAM_MAPPED(SAMTOOLS.out.bambai)
        bambai_ch = FILTER_BAM_MAPPED.out.filtered_bambai
    }
    else {
        bambai_ch = SAMTOOLS.out.bambai
    }

    emit:
    bambai = bambai_ch // sorted and indexed reads (optionally filtered to mapped only): [val(meta), path(bam), path(bai)]
    flagstat = SAMTOOLS.out.flagstat // alignment flagstats [val(meta), path(flagstat_file)]
    cramino_stats = CRAMINO.out.stats_txt // cramino bam stats [val(meta), path(cramino_stat_file)]
    nanostats = NANOPLOT.out.txt // nanostats [val(meta), path(nanostat_file)]
    rseqc_read_dist = rseqc_read_dist_ch // RSeQC read distribution calculations [val(meta), path(read_dist_file)] or Channel.empty(), if skip_rseqc
}
