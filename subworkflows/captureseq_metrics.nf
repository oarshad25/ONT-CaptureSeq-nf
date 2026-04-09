#!/usr/bin/env nextflow

/*
This subworkflow computes metrics for targeted nanopore CaptureSeq data including on-target/off-target counts for loci
in the capture panel and read stats (count, N50, transcript-boundary-spanning count) for genes of particular interest and produces outputs to be
used by main workflow multiqc process.

- Input requirements:
    * BAMs must be coordinate-sorted, indexed
    * panel_gene_list and goi_gene_list are plain text files with one gene
      ID or name per line.
        * Duplicates are silently discarded.
        * GOI gene list should be a subset of the panel gene list.

- Metrics produced:
    * Panel-level:
        * On/off-target read count based on overlap with panel gene-body intervals (bargraph, MultiQC custom content):
          Indicates capture efficiency.
    * GOI-level:
        * Read count per GOI locus and Read N50 per GOI locus: Read count indicates the amount of sequence support
          observed over a GOI locus; N50 indicates whether those reads are long enough to support isoform-level interpretation.
          Together these are practical indicators of whether a GOI has sufficient data for misoform discovery.
        * Count of reads spanning the annotated start and end of the longest transcript per GOI:
          Computes the proportion of GOI reads that span the full length of the longest annotated transcript for each gene.
          This is the most direct indicator of whether Nanopore CaptureSeq is achieving full-length isoform capture at the GOIs,
          which is the primary advantage of the technology over short-read approaches.

- Outputs:
    * panel.bed / goi.bed: BED files derived from input gene lists
    * panel_not_found.txt / goi_not_found.txt: genes absent from GTF (optional, only emitted if missing genes exist)
    * multiqc_files: collected QC files for consumption by the main workflow MultiQC process

- Outline:
    * Create a BED file of the genomic coordinates of genes in the capture panel
    * Calculate number of reads that overlap capture-panel gene-body intervals (on-target) using samtools against panel BED and output
      a tsv to be rendered as a plot in multiqc showing the on/off-target split per sample
    * Create a BED file of the genomic coordinates of genes of particular interest (GOI) if a list of such genes is provided as input
    * Calculate the following metrics per GOI:
        * read count and read N50: The number of reads and their N50 length that overlap each GOI locus
        * transcript-boundary-spanning read count for GOIs. These reads are defined as those whose aligned reference span covers
          the annotated start and end of the longest transcript for a gene, allowing for some tolerance defined by 'goi_fulllength_margin'
          so reads can miss the transcript start and end by a certain number of reference bases and still be counted.
*/

workflow CAPTURESEQ_METRICS {
    take:
    bambai // channel, sorted and indexed bams: [val(meta), path(bam), path(bai)]
    annotation // path, genome annotation as GTF
    panel_gene_list // path, plain text file with list of genes in the capture panel, one geneID/name per line
    goi_gene_list // path, optional plain text file with few select genes of particular interest whose read stats are to be included in final multiQC, one geneID/name per line
    goi_fulllength_margin // int, margin parameter for counting GOI reads as spanning transcript boundaries — number of aligned reference bases that a read can miss at the annotated transcript start and end and still be counted

    main:

    // initialise channel to collect multiqc files generated in this subworkflow
    multiqc_files_ch = Channel.empty()

    // create BED file of coordinates of genes in target capture panel
    MAKE_PANEL_BED(annotation, panel_gene_list)

    // Surface partially missing genes as a pipeline warning because
    // stderr from successful tasks is not reliably shown in the console.
    MAKE_PANEL_BED.out.not_found.subscribe { not_found_file ->
        log.warn("Some genes in ${panel_gene_list.name} were not found in the GTF. See ${not_found_file} for the genes in the capture panel absent in the GTF.")
    }

    // Calculate number of primary mapped reads that fall within the capture panel gene list.
    ON_TARGET_COUNTS(bambai, MAKE_PANEL_BED.out.bed)

    multiqc_files_ch = multiqc_files_ch.mix(ON_TARGET_COUNTS.out.mqc_tsv.map { _meta, f -> f })

    // summarise per-sample QC metrics for reads overlapping each GOI locus if a GOI list is provided
    if (goi_gene_list.name != 'NO_FILE') {
        // Build a BED of the curated genes whose read support metrics should be highlighted in MultiQC.
        MAKE_GOI_BED(annotation, goi_gene_list)

        // Surface partially missing genes as a pipeline warning because
        // stderr from successful tasks is not reliably shown in the console.
        MAKE_GOI_BED.out.not_found.subscribe { not_found_file ->
            log.warn("Some genes in ${goi_gene_list.name} were not found in the GTF. See ${not_found_file} for the genes of interest absent in the GTF.")
        }

        // calculate read counts and read N50 for genes of interest
        GOI_READCOUNT_N50(bambai, MAKE_GOI_BED.out.bed)

        // calculate the number of reads whose aligned reference span covers the annotated start and end
        // of the longest transcript for each gene of interest, within the configured margin
        GOI_FULLLENGTH_COUNT(bambai, MAKE_GOI_BED.out.bed, annotation, goi_fulllength_margin)

        // Join per-GOI metric channels by sample meta before combining into
        // a single MultiQC TSV — ensures metrics are correctly paired per sample
        goi_tsvs_ch = GOI_READCOUNT_N50.out.tsv.join(GOI_FULLLENGTH_COUNT.out.tsv, by: 0)

        GOI_COMBINE_MQC(goi_tsvs_ch)

        multiqc_files_ch = multiqc_files_ch.mix(GOI_COMBINE_MQC.out.mqc_tsv.map { _meta, f -> f })
    }

    emit:
    panel_bed = MAKE_PANEL_BED.out.bed // panel.bed BED file of genomic coordinates of genes in the capture panel
    panel_not_found = MAKE_PANEL_BED.out.not_found.ifEmpty([]) // panel_not_found.txt or empty channel, list of genes in the capture panel that were not found in the GTF
    goi_bed = MAKE_GOI_BED.out.bed //goi.bed, BED file of genomic coordinates of genes of interest
    goi_not_found = MAKE_GOI_BED.out.not_found.ifEmpty([]) // goi_not_found.txt or empty channel, list of genes of interest that were not found in the GTF
    multiqc_files = multiqc_files_ch // all QC files for main workflow MultiQC input, including panel locus-overlap counts and GOI metrics if a GOI list is provided
}


// Converts the capture panel gene list to a BED file of gene body
// coordinates by looking up each gene in the GTF. Ensembl version
// suffixes are stripped from GTF gene_id fields before matching.
// Genes not found in the GTF are written to panel_not_found.txt;
// this file is only emitted if at least one gene is missing.
process MAKE_PANEL_BED {
    label 'low'

    conda 'conda-forge::python=3.11'
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/python:3.11'
        : 'quay.io/biocontainers/python:3.11'}"

    input:
    path gtf
    path panel_gene_list

    output:
    path "panel.bed", emit: bed
    path "panel_not_found.txt", optional: true, emit: not_found

    script:
    """
    make_genes_bed.py \\
        --gtf       ${gtf} \\
        --gene-list ${panel_gene_list} \\
        --output    panel.bed \\
        --not-found panel_not_found.txt
    """
}

// Converts the genes of interest gene list to a BED file of gene body coordinates.
// Identical logic to MAKE_PANEL_BED but kept as a separate process so
// either can be modified, retried, or removed independently.
process MAKE_GOI_BED {
    label 'low'

    conda 'conda-forge::python=3.11'
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/python:3.11'
        : 'quay.io/biocontainers/python:3.11'}"

    input:
    path gtf
    path goi_gene_list

    output:
    path "goi_genes.bed", emit: bed
    path "goi_genes_not_found.txt", optional: true, emit: not_found

    script:
    """
    make_genes_bed.py \\
        --gtf       ${gtf} \\
        --gene-list ${goi_gene_list} \\
        --output    goi_genes.bed \\
        --not-found goi_genes_not_found.txt
    """
}

// Computes on-target and off-target read counts for each sample using
// samtools view -c against the panel BED. Outputs a MultiQC bargraph
// custom content file showing the on/off-target split per sample.
// The cpswitch option allows toggling between raw counts and percentage
// views in the MultiQC report.
process ON_TARGET_COUNTS {
    tag "${meta.id}"
    label 'low'

    conda 'bioconda::samtools=1.21'
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.21--h96c455f_1'
        : 'quay.io/biocontainers/samtools:1.21--h96c455f_1'}"

    input:
    tuple val(meta), path(bam), path(bai)
    path panel_bed

    output:
    tuple val(meta), path("*_ontarget_mqc.tsv"), emit: mqc_tsv

    script:
    """
    # total primary mapped reads
    TOTAL=\$(samtools view -c -@ ${task.cpus} -F 2308 ${bam})

    # primary mapped reads that overlap panel gene-body BED intervals
    ON=\$(samtools view -c -@ ${task.cpus} -F 2308 -L ${panel_bed} ${bam})

    # off-target reads are the remainder of the total primary mapped reads
    OFF=\$((TOTAL - ON))

    cat > ${meta.id}_ontarget_mqc.tsv << EOF
# id: 'on_target'
# section_name: 'Panel On-target'
# description: 'Primary mapped reads overlapping capture-panel gene-body intervals'
# format: 'tsv'
# plot_type: 'bargraph'
# pconfig:
#   id: 'on_target_counts_plot'
#   title: 'On-target'
#   ylab: 'Reads'
#   cpswitch: true
#   cpswitch_c_active: false
# categories:
#   on_target:
#     name: 'On-target Reads'
#     color: '#2f7ed8'
#   off_target:
#     name: 'Off-target Reads'
#     color: '#d9534f'
Sample\ton_target\toff_target
${meta.id}\t\${ON}\t\${OFF}
EOF
    """
}

// Counts primary mapped reads and their N50 read lengthper gene of interest.
process GOI_READCOUNT_N50 {
    tag "${meta.id}"
    label 'low'

    conda 'bioconda::pysam=0.22'
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.22.1--py39hdd5828d_3'
        : 'quay.io/biocontainers/pysam:0.22.1--py39hdd5828d_3'}"

    input:
    tuple val(meta), path(bam), path(bai)
    path goi_bed

    output:
    tuple val(meta), path("*_goi_readcount_n50.tsv"), emit: tsv

    script:
    """
    goi_readcount_n50.py \\
        --bam    ${bam} \\
        --goi    ${goi_bed} \\
        --sample ${meta.id} \\
        --output ${meta.id}_goi_readcount_n50.tsv
    """
}

// Computes the number of GOI reads whose aligned reference span covers the
// annotated start and end of the longest transcript for each gene.
process GOI_FULLLENGTH_COUNT {
    tag "${meta.id}"
    label 'low'

    conda 'bioconda::pysam=0.22'
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.22.1--py39hdd5828d_3'
        : 'quay.io/biocontainers/pysam:0.22.1--py39hdd5828d_3'}"

    input:
    tuple val(meta), path(bam), path(bai)
    path goi_bed
    path gtf
    val goi_fulllength_margin

    output:
    tuple val(meta), path("*_goi_fulllength.tsv"), emit: tsv

    script:
    """
    goi_fulllength_count.py \\
        --bam    ${bam} \\
        --goi    ${goi_bed} \\
        --gtf    ${gtf} \\
        --margin ${goi_fulllength_margin} \\
        --sample ${meta.id} \\
        --output ${meta.id}_goi_fulllength.tsv
    """
}


// Merges the per-GOI metric TSVs from GOI_READCOUNT_N50,
// GOI_FULLLENGTH_COUNT into a single MultiQC table section per sample
process GOI_COMBINE_MQC {
    tag "${meta.id}"
    label 'single'

    conda 'conda-forge::python=3.11'
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/python:3.11'
        : 'quay.io/biocontainers/python:3.11'}"

    input:
    tuple val(meta), path(readcount_n50_tsv), path(fulllength_tsv)

    output:
    tuple val(meta), path("*_goi_summary_mqc.tsv"), emit: mqc_tsv

    script:
    """
    goi_combine_mqc.py \\
        --readcount-n50 ${readcount_n50_tsv} \\
        --fulllength    ${fulllength_tsv} \\
        --sample        ${meta.id} \\
        --output        ${meta.id}_goi_summary_mqc.tsv
    """
}
