#!/usr/bin/env nextflow

// Preprocess reference files

// This subworkflow takes input reference fasta and GTF file paths and decompresses
// them if needed. It converts the input GTF to BED as well using the UCSC utilities
// gtfToGenePred and genePredToBed for downstream subworkflows/processes that need BED input.
// Output is reference file channels of uncompressed fasta and annotation as GTF and BED


// TODO: Replace gzip with pigz for fast decompression
// to unzip annotation if needed
include { GUNZIP as GUNZIP_GENOME } from '../modules/gunzip'
include { GUNZIP as GUNZIP_ANNOTATION } from '../modules/gunzip'

// these UCSC utilities are used to convert the input GTF to BED
// The output BED is then fed to RSEQC subworkflow as it takes annotation in BED format
include { UCSC_GTF_TO_GENEPRED } from '../modules/ucsc_gtf_to_genepred.nf'
include { UCSC_GENEPRED_TO_BED } from '../modules/ucsc_genepred_to_bed.nf'

workflow PREPARE_REFERENCE_FILES {
    take:
    genome // path/to/genome/fasta
    annotation // path/to/annotation/gtf

    main:

    // unzip genome fasta if compressed
    if (file(genome).extension == "gz") {
        // unzip genome fasta
        GUNZIP_GENOME(genome)
        genome_ch = GUNZIP_GENOME.out.gunzip_file
    }
    else {
        genome_ch = file(genome)
    }

    // unzip genome gtf if compressed
    if (file(annotation).extension == "gz") {
        GUNZIP_ANNOTATION(annotation)
        annotation_ch = GUNZIP_ANNOTATION.out.gunzip_file
    }
    else {
        annotation_ch = annotation
    }

    // convert GTF to BED using UCSC utilities:
    // * First convert the input GTF to genePred by using UCSC's
    //   gtfToGenePred utility.
    // * The genePred is then converted to a BED using UCSC's genePredToBed.

    UCSC_GTF_TO_GENEPRED(annotation_ch)
    genepred_ch = UCSC_GTF_TO_GENEPRED.out.genepred

    UCSC_GENEPRED_TO_BED(genepred_ch)
    annotation_bed_ch = UCSC_GENEPRED_TO_BED.out.bed

    emit:
    genome_ch = genome_ch // path to genome fasta
    annotation_ch = annotation_ch // path to gtf
    annotation_bed_ch = annotation_bed_ch // path to annotation as bed
}
