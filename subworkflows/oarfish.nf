/*
Oarfish transcript quantification subworkflow

Builds one transcriptome index, quantifies each sample's QC'ed reads, and
merges per-sample .quant files into a transcript x sample counts matrix.
*/

workflow OARFISH {
    take:
    reads // [ val(meta), path(fastq) ]
    transcriptome // path to transcriptome fasta
    seq_tech // Oarfish --seq-tech option
    filter_group // Oarfish --filter-group option
    model_coverage // Boolean, whether to pass --model-coverage
    extra_opts // string, any additional options to pass to Oarfish quantification

    main:
    OARFISH_INDEX(transcriptome, seq_tech)

    OARFISH_QUANT(
        reads,
        OARFISH_INDEX.out.index,
        seq_tech,
        filter_group,
        model_coverage,
        extra_opts,
    )

    OARFISH_MERGE_QUANTS(
        OARFISH_QUANT.out.quant.map { _meta, quant -> quant }.collect()
    )

    emit:
    index = OARFISH_INDEX.out.index
    quant = OARFISH_QUANT.out.quant
    meta_info = OARFISH_QUANT.out.meta_info
    ambig_info = OARFISH_QUANT.out.ambig_info
    infreps = OARFISH_QUANT.out.infreps
    assignment_probs = OARFISH_QUANT.out.assignment_probs
    assignment_probs_lz4 = OARFISH_QUANT.out.assignment_probs_lz4
    transcript_counts_matrix = OARFISH_MERGE_QUANTS.out.transcript_counts_matrix
}

process OARFISH_INDEX {
    label "medium"

    publishDir "${params.outdir}/oarfish/index", mode: "link"

    conda "bioconda::oarfish=0.9.4--h7f5d12c_0"
    container "quay.io/biocontainers/oarfish:0.9.4--h7f5d12c_0"

    input:
    // transcriptome fasta
    path transcriptome

    // Oarfish --seq-tech option
    val seq_tech

    output:
    path "transcriptome.mmi", emit: index

    script:
    """
    oarfish \\
        --only-index \\
        --annotated ${transcriptome} \\
        --seq-tech ${seq_tech} \\
        --threads ${task.cpus} \\
        --index-out transcriptome.mmi
    """

    stub:
    """
    touch transcriptome.mmi
    """
}

process OARFISH_QUANT {
    tag "${meta.id}"

    label "medium"

    publishDir "${params.outdir}/oarfish/${meta.id}", mode: "link"

    conda "bioconda::oarfish=0.9.4--h7f5d12c_0"
    container "quay.io/biocontainers/oarfish:0.9.4--h7f5d12c_0"

    input:
    // QC'ed reads
    tuple val(meta), path(reads)

    // prebuilt transcriptome index
    path index

    // Oarfish --seq-tech option
    val seq_tech

    // Oarfish --filter-group option
    val filter_group

    // Boolean, whether to pass --model-coverage
    val model_coverage

    // any additional command line options
    val extra_opts

    output:
    tuple val(meta), path("${meta.id}.quant"), emit: quant
    tuple val(meta), path("${meta.id}.meta_info.json"), emit: meta_info
    tuple val(meta), path("${meta.id}.ambig_info.tsv"), emit: ambig_info
    tuple val(meta), path("${meta.id}.infreps.pq"), emit: infreps, optional: true
    tuple val(meta), path("${meta.id}.prob"), emit: assignment_probs, optional: true
    tuple val(meta), path("${meta.id}.prob.lz4"), emit: assignment_probs_lz4, optional: true

    script:
    def model_coverage_flag = model_coverage ? "--model-coverage" : ""

    """
    oarfish \\
        --reads ${reads} \\
        --index ${index} \\
        --seq-tech ${seq_tech} \\
        --threads ${task.cpus} \\
        --output ${meta.id} \\
        --filter-group ${filter_group} \\
        ${model_coverage_flag} \\
        ${extra_opts}
    """

    stub:
    """
    printf "target_id\\tnum_reads\\ntranscript_1\\t0\\n" > ${meta.id}.quant
    printf "{}\\n" > ${meta.id}.meta_info.json
    printf "target_id\\tunique\\tambiguous\\ttotal\\ntranscript_1\\t0\\t0\\t0\\n" > ${meta.id}.ambig_info.tsv
    """
}

process OARFISH_MERGE_QUANTS {
    label "low"

    publishDir "${params.outdir}/oarfish/", mode: "copy"

    conda "conda-forge::python=3.11"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/python:3.11'
        : 'quay.io/biocontainers/python:3.11'}"

    input:
    // list of Oarfish .quant files
    path quant_files

    output:
    path "oarfish.transcript_counts.tsv", emit: transcript_counts_matrix

    script:
    """
    oarfish_merge_quants.py \\
        --output oarfish.transcript_counts.tsv \\
        ${quant_files}
    """

    stub:
    """
    printf "transcript_id\\n" > oarfish.transcript_counts.tsv
    """
}
