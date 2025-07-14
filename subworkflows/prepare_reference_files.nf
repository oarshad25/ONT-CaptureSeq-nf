#!/usr/bin/env nextflow

// Preprocess reference files

// Create reference file channels from input strings of paths to reference files.
// * Uncompress reference files if they are gzipped
// Output is path of processed references

// TODO: Replace gzip with pigz for fast decompression

include { GUNZIP as GUNZIP_GENOME } from '../modules/gunzip'
include { GUNZIP as GUNZIP_ANNOTATION } from '../modules/gunzip'

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
        annotation_ch = file(annotation)
    }

    emit:
    genome_ch = genome_ch // path to genome fasta
    annotation_ch = annotation_ch // path to gtf
}
