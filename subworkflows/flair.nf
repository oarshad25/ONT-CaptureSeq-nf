/*
FLAIR subworkflow

Subworkflow to run FLAIR to generate a transcriptome of high confidence
isoforms from aligned bams

*/

include { MERGE_BAMS } from "../modules/merge_bams"

workflow FLAIR {
    take:
    merged_bambai // merged and indexed bam: tuple path(bam), path(bai)
    genome // path to genome fasta
    annotation // path to annotation gtf
    reads // path to fastq reads, list of fastq files
    flair_collapse_extra_opts // string, any additional options to pass to flair_collapse module
    flair_align_reads_manifest // path to reads manifest for flair align module

    main:
    FLAIR_BAM2BED12(merged_bambai)

    FLAIR_CORRECT(
        FLAIR_BAM2BED12.out.bed,
        genome,
        annotation,
    )

    FLAIR_COLLAPSE(
        FLAIR_CORRECT.out.corrected_bed,
        genome,
        annotation,
        reads,
        flair_collapse_extra_opts,
    )

    FLAIR_QUANTIFY(
        FLAIR_COLLAPSE.out.isoforms_fa,
        flair_align_reads_manifest,
        reads,
    )

    emit:
    isoforms_gtf = FLAIR_COLLAPSE.out.isoforms_gtf // path, isoforms gtf file
    isoforms_bed = FLAIR_COLLAPSE.out.isoforms_bed // path, isoforms bed file
    isoforms_fa = FLAIR_COLLAPSE.out.isoforms_fa // path, isoforms fasta file
    transcript_counts_tsv = FLAIR_QUANTIFY.out.flair_quantify_tsv // path, isoform quantification tsv
}

// convert merged bam to bed12
process FLAIR_BAM2BED12 {
    label "low"

    conda "bioconda::flair=2.2.0"
    container "quay.io/biocontainers/flair:2.2.0--pyhdfd78af_0"

    input:
    // path to merged bam
    tuple path(bam), path(bai)

    output:
    path "merged.bed", emit: bed

    script:
    """
    bam2Bed12 -i ${bam} > merged.bed
    """
}

// run flair correct module on merged bed to correct misaligned splice sites using genome annotation
process FLAIR_CORRECT {
    label "medium"

    conda "bioconda::flair=2.2.0"
    container "quay.io/biocontainers/flair:2.2.0--pyhdfd78af_0"

    input:
    // path to merged bed
    path bed
    // genome fasta file
    path genome
    // annotation gtf
    path annotation

    output:
    path "flair_all_corrected.bed", emit: corrected_bed
    path "flair_all_inconsistent.bed", emit: inconsistent_bed, optional: true
    path "flair_cannot_verify.bed", emit: cannot_verify_bed, optional: true

    script:
    """
    flair correct \\
        --threads ${task.cpus} \\
        --query ${bed} \\
        --genome ${genome} \\
        --gtf ${annotation}
    """
}

// run flair collapse module to define high-confidence isoforms from corrected reads
process FLAIR_COLLAPSE {
    label "high"

    publishDir "${params.outdir}/flair", mode: 'link'

    conda "bioconda::flair=2.2.0"
    container "quay.io/biocontainers/flair:2.2.0--pyhdfd78af_0"

    input:
    // path to corrected bed output from flair_correct
    path bed
    // genome fasta file
    path genome
    // annotation gtf
    path annotation
    // list of raw read fastqs
    path reads
    // any additional command line options
    val extra_opts

    output:
    path "flair_collapse", emit: flair_collapse_out_dir

    path ("**/flair.collapse.isoforms.fa"), emit: isoforms_fa
    path ("**/flair.collapse.isoforms.bed"), emit: isoforms_bed
    path ("**/flair.collapse.isoforms.gtf"), emit: isoforms_gtf
    path ("**/flair.collapse.isoform.read.map.txt"), emit: isoforms_read_map_txt, optional: true

    script:
    """
    mkdir -p flair_collapse

    flair collapse \\
        --threads ${task.cpus} \\
        --query ${bed} \\
        --genome ${genome} \\
        --gtf ${annotation} \\
        --reads ${reads} \\
        --output flair_collapse/flair.collapse \\
        ${extra_opts}
    """
}

process FLAIR_QUANTIFY {
    label "medium"

    publishDir "${params.outdir}/flair/flair_quantify", mode: 'link'

    conda "bioconda::flair=2.2.0"
    container "quay.io/biocontainers/flair:2.2.0--pyhdfd78af_0"

    input:
    // path to isoforms fasta output from flair collapse
    path isoforms
    // path to reads manifest
    path reads_manifest
    // path to raw read fastqs
    path reads

    output:
    path "flair.quantify.counts.tsv", emit: flair_quantify_tsv

    script:
    """
    # rename fastq files to just have meta.id by removing anything after a . to match read manifest
    for fq in *.fastq
    do
        mv "\$fq" "\${fq%%.*}.fastq"
    done

    flair quantify \\
        --threads ${task.cpus} \\
        --isoforms ${isoforms} \\
        --reads_manifest ${reads_manifest} \\
        --sample_id_only
    """
}
