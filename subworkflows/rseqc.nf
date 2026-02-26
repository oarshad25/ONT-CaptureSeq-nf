#!/usr/bin/env nextflow

/*
 * Alignment QC with RSeQC

Use RSeQC on aligned reads. This subworkflow takes the aligned reads
(sorted and indexed BAMs) and input reference gene model bed and runs the RSeqC modules
read_gc.py, read_distribution.py, and infer_experiment.py
to calculate various statistics

* RSeQC script read_gc.py is used to calculate GC content distribution of reads
* RSeQC module read_distribution.py is used to calculate read distribution
  of mapped reads over genomic features.
  to perform alternative splicing analysis
* Module infer_experiment.py is used to infer strandedness, guess how RNA-seq experiment
  was configured
*/
workflow RSEQC {
    take:
    bambai // sorted and indexed aligned reads: [val(meta), path(bam), path(bai)]
    rseqc_bed_ch // path/to/bed: path(bed)

    main:
    RSEQC_READ_GC(bambai)
    RSEQC_READ_DISTRIBUTION(bambai, rseqc_bed_ch)
    RSEQC_INFER_EXPERIMENT(bambai, rseqc_bed_ch)

    emit:
    read_gc_xls = RSEQC_READ_GC.out.xls // [val(meta), path(xls)]
    read_distribution_txt = RSEQC_READ_DISTRIBUTION.out.txt // [val(meta), path(txt)]
    infer_experiment_txt = RSEQC_INFER_EXPERIMENT.out.txt // [val(meta), path(txt)]
}

// use RSeQC to compute GC content distribution of reads
process RSEQC_READ_GC {
    tag "${meta.id}"
    label 'single'

    conda "bioconda rseqc=5.0.4"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/rseqc:5.0.4--pyhdfd78af_0'
        : 'quay.io/biocontainers/rseqc:5.0.4--pyhdfd78af_0'}"

    input:
    // sorted and indexed bam file
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.xls"), emit: xls
    tuple val(meta), path("*.r"), emit: rscript, optional: true
    tuple val(meta), path("*.pdf"), emit: pdf, optional: true

    script:
    """
    read_GC.py \\
        -i ${bam} \\
        -o ${meta.id}
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

// use RSeQC to infer strandedness
process RSEQC_INFER_EXPERIMENT {
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
    tuple val(meta), path("${meta.id}.infer_experiment.txt"), emit: txt

    script:
    """
    infer_experiment.py \\
        -i ${bam} \\
        -r ${bed} \\
        > ${meta.id}.infer_experiment.txt
    """
}
