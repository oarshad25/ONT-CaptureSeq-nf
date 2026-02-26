#!/usr/bin/env nextflow

/*
 * Alignment QC with RSeQC

Use RSeQC on aligned reads. This subworkflow takes the aligned reads
(sorted and indexed BAMs) and input reference gene model bed and runs the RSeqC modules
read_gc.py, read_distribution.py, infer_experiment.py and geneBody_coverage.py
to calculate various statistics

* RSeQC script read_gc.py is used to calculate GC content distribution of reads
* RSeQC module read_distribution.py is used to calculate read distribution
  of mapped reads over genomic features.
  to perform alternative splicing analysis
* Module infer_experiment.py is used to infer strandedness, guess how RNA-seq experiment
  was configured
* Module geneBody_coverage.py is used to calculate gene body coverage of mapped reads.
  This module is only run if a bed file containing housekeeping genes is provided (parameter 'rseqc_housekeeping_bed')
*/
workflow RSEQC {
    take:
    bambai // sorted and indexed aligned reads: [val(meta), path(bam), path(bai)]
    rseqc_bed_ch // path/to/bed: path(bed), annotation in bed format
    rseqc_housekeeping_bed_ch // path/to/bed: path(bed), housekeeping genes bed for gene body coverage module

    main:
    RSEQC_READ_GC(bambai)
    RSEQC_READ_DISTRIBUTION(bambai, rseqc_bed_ch)
    RSEQC_INFER_EXPERIMENT(bambai, rseqc_bed_ch)

    // run RSeQC gene body coverage module if housekeeping gene bed file is provided
    rseqc_genebody_coverage_txt_ch = Channel.empty()
    if (rseqc_housekeeping_bed_ch.name != 'NO_FILE') {
        RSEQC_GENEBODY_COVERAGE(bambai, rseqc_housekeeping_bed_ch)
        rseqc_genebody_coverage_txt_ch = RSEQC_GENEBODY_COVERAGE.out.txt
    }

    emit:
    read_gc_xls = RSEQC_READ_GC.out.xls // [val(meta), path(xls)]
    read_distribution_txt = RSEQC_READ_DISTRIBUTION.out.txt // [val(meta), path(txt)]
    infer_experiment_txt = RSEQC_INFER_EXPERIMENT.out.txt // [val(meta), path(txt)]
    rseqc_genebody_coverage_txt = rseqc_genebody_coverage_txt_ch
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

// use RSeQC to compute gene body coverage
process RSEQC_GENEBODY_COVERAGE {
    tag "${meta.id}"
    label 'low'

    conda "bioconda rseqc=5.0.4"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/rseqc:5.0.4--pyhdfd78af_0'
        : 'quay.io/biocontainers/rseqc:5.0.4--pyhdfd78af_0'}"

    input:
    // sorted and indexed bam file
    tuple val(meta), path(bam), path(bai)

    // housekeeping genes bed file for gene body coverage calculation
    path bed

    output:
    tuple val(meta), path("${meta.id}.geneBodyCoverage.txt"), emit: txt

    script:
    """
    geneBody_coverage.py \\
        -r ${bed} \\
        -i ${bam} \\
        -o ${meta.id}
    """
}
