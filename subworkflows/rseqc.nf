#!/usr/bin/env nextflow

/*
 * Alignment QC with RSeQC

Use RSeQC on aligned reads. This subworkflow takes the aligned reads
(sorted and indexed BAMs) and input reference GTF and calculates how mapped
reads are distributed over genome features and outputs text files containing
read distribution calculations

* First it converts the input GTF to genePred by using UCSC's
  gtfToGenePred utility.
* The genePred is then converted to a BED using UCSC's genePredToBed.
* The reference gene model BED file along with the mapped reads are passed as
  input to the RSeQC script read_distribution.py to calculate read distribution
  of mapped reads over genomic features.
*/
workflow RSEQC {
    take:
    bambai // sorted and indexed aligned reads: [val(meta), path(bam), path(bai)]
    gtf // path/to/gtf: path(gtf)

    main:
    UCSC_GTF_TO_GENEPRED(gtf)
    genepred_ch = UCSC_GTF_TO_GENEPRED.out.genepred

    UCSC_GENEPRED_TO_BED(genepred_ch)
    rseqc_bed_ch = UCSC_GENEPRED_TO_BED.out.bed

    RSEQC_READ_DISTRIBUTION(bambai, rseqc_bed_ch)

    emit:
    txt = RSEQC_READ_DISTRIBUTION.out.txt // [val(meta), path(txt)]
}

// use UCSC's gtfToGenePred to convert a GTF file to genePred
process UCSC_GTF_TO_GENEPRED {
    label "single"

    conda "bioconda ucsc-gtftogenepred=469"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/ucsc-gtftogenepred:469--h664eb37_1'
        : 'quay.io/biocontainers/ucsc-gtftogenepred:469--h664eb37_1'}"

    input:
    path gtf

    output:
    path "*.genepred", emit: genepred

    script:
    """
    gtfToGenePred \\
        ${gtf}  \\
        ${gtf.simpleName}.genepred
    """
}

// use UCSC's genePredToBed to convert genepred to bed
process UCSC_GENEPRED_TO_BED {
    label "single"

    conda "bioconda ucsc-genepredtobed=469"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/ucsc-genepredtobed:469--h664eb37_1'
        : 'quay.io/biocontainers/ucsc-genepredtobed:469--h664eb37_1'}"

    input:
    path genepred

    output:
    path "*.bed", emit: bed

    script:
    """
    genePredToBed \\
        ${genepred}  \\
        ${genepred.simpleName}.bed
    """
}

// use RSeQC to calculate read distribution of aligned reads over genomic features
process RSEQC_READ_DISTRIBUTION {
    tag "${meta.id}"
    label 'single'

    conda "bioconda rseqc=5.0.4"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/rseqc:5.0.4--pyhdfd78af_0'
        : 'quay.io/biocontainers/rseqc:5.0.4--pyhdfd78af_0'}"

    input:
    // sorted and indexed bam file
    tuple val(meta), path(bam), path(bai)

    // bed file output from UCSC utilities
    path bed

    output:
    tuple val(meta), path("${meta.id}.read_distribution.txt"), emit: txt

    script:
    """
    read_distribution.py \\
        -i ${bam} \\
        -r ${bed} \\
        > ${meta.id}.read_distribution.txt
    """
}
