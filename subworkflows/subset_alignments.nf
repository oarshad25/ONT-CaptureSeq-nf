#!/usr/bin/env nextflow

/*
Subset alignments subworkflow

* Subset filtered alignments to genes of interest:
* Takes input annotation GTF and genelist and outputs indexed bams subset to genes in list
* 1. First creates a BED file from the GTF containing the genomic coordinates of genes of
* interest in BED format from GTF.
* 2. Then uses samtools to subset the input aligned reads to the
* coordinates in the BED file

*/

workflow SUBSET_ALIGNMENTS {
    take:
    bambai // sorted and indexed bam: [val(meta), path(bam), path(bai)]
    annotation // annotation as GTF
    genelist // file with ENSEMBL ID's of genes of interest to filter input bams to

    main:

    // take the annotation GTF and genelist and create a BED file
    // with the genomic coordinates of genes of interest
    CREATE_GENES_OF_INTEREST_BED(annotation, genelist)

    bed_ch = CREATE_GENES_OF_INTEREST_BED.out.bed

    // subset input alignment bam to regions in BED
    SUBSET_BAM_WITH_BED(bambai, bed_ch)

    // reads aligned to genes of interest
    bambai_subset_ch = SUBSET_BAM_WITH_BED.out.bambai_subset

    emit:
    bambai = bambai_subset_ch // reads aligned to genes of interest: [val(meta), path(bam), path(bai)]
}

// take the annotation GTF and genelist and create a BED file
// with the genomic coordinates of genes of interest
process CREATE_GENES_OF_INTEREST_BED {
    label "single"
    label "local_software"

    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:24.04'
        : 'quay.io/biocontainers/ubuntu:24.04'}"

    input:
    // path to annotation GTF
    path annotation
    // path to genelist to filter GTF to only these genes
    path genelist

    output:
    // GTF subset to genes of interest
    path "genes_of_interest.bed", emit: bed

    script:
    """
    # 1. Subset GTF
    # subset GTF to genes of interest
    # -f option is used to search on file
    # -F for eact matches

    grep -Ff ${genelist} ${annotation} > genes_of_interest.gtf

    # 2. Convert to BED
    # extract the gene entries and convert to four column BED with chromosome, 0 based start index, end, gene_id
    # filters the GTF for gene entries and extracts 1st column (chromosome),
    # subtracts one from the fourth column which is start position (GTF uses 1 based index)
    # takes 5th column which is end position
    # since awk seperates columns based on space, gene ID is 10th column e.g. "ENSG00000012660";
    # use tr -d to delete " and ; characters from ID

    awk '\$3 == "gene" {print \$1"\t"\$4-1"\t"\$5"\t"\$10}' genes_of_interest.gtf | \
        tr -d '";' > genes_of_interest.bed
    """
}

// use samtools to filter sorted and indexed aligned BAM to regions specified in input BED
process SUBSET_BAM_WITH_BED {
    tag "${meta.id}"
    label 'low'

    publishDir "${params.outdir}/bam/genes_of_interest/", mode: "link"

    conda "bioconda samtools=1.21"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.21--h96c455f_1'
        : 'quay.io/biocontainers/samtools:1.21--h96c455f_1'}"

    input:
    // sorted and indexed bam file
    tuple val(meta), path(bam), path(bai)

    // bed file giving regions to filter the input bams to
    path bed

    output:
    // filtered and indexed bam file overlapping regions specified in input bed
    tuple val(meta), path("${meta.id}.subset.bam"), path("${meta.id}.subset.bam.bai"), emit: bambai_subset

    script:
    """
    # filter input bam to regions given in bed
    samtools view -@ ${task.cpus} -bh -L ${bed} ${bam} > ${meta.id}.subset.bam

    # index filtered BAM
    samtools index -@ ${task.cpus} ${meta.id}.subset.bam
    """
}
