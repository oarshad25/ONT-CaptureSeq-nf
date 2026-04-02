#!/usr/bin/env nextflow

/*

- This subworkflow computes:

* on-target rate (Capture Efficiency): Fraction of reads mapping to genes in target panel
* Stats for genes of particular interest: Read stats (total and N50) for a few select genes of partcular interest specified in a list if such a list is provided

and produces tables to be included in multiqc

- Outline:
* Calculate on-target metrics using the gene counts matrix and list of panel genes using a custom script outputting a table for multiQC
* Create a BED file from the GTF containing the genomic coordinates of select genes of interest specified in the optional input if provided using GNU awk
  If any genes in teh select gene list are not present in the GTF surface this as a warning and point to teh file listing the missing genes
  Error out altogther if none of the genes listed in the select gene list are present in the GTF.
* Calculate read stats for these genes using a custom script that takes the BED file of select genes produced above and
  input of indexed BAMs producing a table for multiqc

*/

workflow PANEL_METRICS {
    take:
    bambai // sorted and indexed bams: [val(meta), path(bam), path(bai)]
    annotation // annotation as GTF
    panel_gene_list // file with list of genes in the capture panel
    select_gene_list // optional file with few select genes of particular interest whose read stats are to be included in final multiQC
    gene_counts_matrix // gene x sample counts matrix

    main:
    // Summarise the fraction of assigned reads that fall within the capture panel gene list.
    CALC_PANEL_ON_TARGET(gene_counts_matrix, panel_gene_list)

    // add output stats table to a multiqc channel
    multiqc_files_ch = CALC_PANEL_ON_TARGET.out.mqc_tsv

    // summarise per-sample read count and N50 for reads overlapping each gene in the select genes list if provided
    if (select_gene_list.name != 'NO_FILE') {
        // Build a BED of the curated genes whose per-gene read stats should be highlighted in MultiQC.
        CREATE_GENES_BED(annotation, select_gene_list)

        // Surface partially missing genes as a pipeline warning because
        // stderr from successful tasks is not reliably shown in the console.
        CREATE_GENES_BED.out.missing_genes.subscribe { missing_genes_file ->
            log.warn("Some genes in ${select_gene_list.name} were not found in the GTF. See ${missing_genes_file} for the genes in the select gene list absent in the GTF.")
        }

        // calculate read stats for the select genes and add the output table to the multiqc channel
        CALC_PANEL_SELECT_GENES_STATS(
            bambai.collect { it -> it[1] },
            bambai.collect { it -> it[2] },
            CREATE_GENES_BED.out.bed,
        )

        // add custom-content multiqc table to multiqc channel.
        multiqc_files_ch = multiqc_files_ch.mix(CALC_PANEL_SELECT_GENES_STATS.out.mqc_tsv)
    }

    emit:
    multiqc_files = multiqc_files_ch // [ path(panel_metrics_on_target_mqc.tsv), optional: path(select_gene_stats_mqc.tsv)] channel of paths to multiqc input files generated in this subworkflow
}

// Convert an input gene list into a BED file of gene intervals from GTF annotation.
// Each BED record is labelled with gene_name when present, otherwise gene_id.
// Error out if none of the genes in input gene list are found in either gene_id or gene_name fields of the GTF annotation
// If only some of the input genes are found in the GTF, surface this as a warning and point to a file listing the missing genes from the input gene list that were not found in the GTF annotation.
process CREATE_GENES_BED {
    label "single"
    label "local_software"

    publishDir "${params.outdir}/panel_metrics/", mode: "copy"

    // need gnu awk (gawk) for 3-argument match function
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/gawk:5.3.0'
        : 'quay.io/biocontainers/gawk:5.3.0'}"

    input:
    path gtf
    path gene_list

    output:
    path "${gtf.baseName}.genes.bed", emit: bed
    path "missing_genes.txt", optional: true, emit: missing_genes

    script:
    """
    # Normalise the input gene list before matching so whitespace,
    # comments, duplicate entries and Ensembl version suffixes do not affect
    # the comparison against the annotation.
    gawk '
        {
            line = \$0
            gsub(/^[[:space:]]+|[[:space:]]+\$/, "", line)
            if (line != "" && line !~ /^#/) {
                sub(/\\.[0-9]+\$/, "", line)
                print line
            }
        }
    ' "${gene_list}" | sort -u > requested_genes.txt

    gawk -v OFS='\t' '
        BEGIN {
            # Load the requested genes into an awk array so we can do fast
            # membership checks while scanning the GTF.
            while ((getline line < "requested_genes.txt") > 0) {
                if (line != "") {
                    genes[line] = 1
                }
            }
            close("requested_genes.txt")
        }

        \$3 == "gene" {
            # Convert GTF gene records into BED coordinates and retain the gene
            # identifiers needed for matching against the requested list.
            chrom     = \$1
            start     = \$4 - 1    # GTF is 1-based; BED is 0-based
            end       = \$5
            gene_id   = ""
            gene_name = ""

            # Pull the gene_id and gene_name attributes out of the GTF column 9
            # annotation field when they are present.
            if (match(\$0, /gene_id "([^"]+)"/, a))   gene_id   = a[1]
            if (match(\$0, /gene_name "([^"]+)"/, b)) gene_name = b[1]

            # strip ENSEMBL version number (e.g., .17)
            sub(/\\.[0-9]+\$/, "", gene_id)

            # Accept either a gene_id or gene_name hit, and prefer gene_name as
            # the BED label when available because it is usually more readable.
            if ((gene_id in genes) || (gene_name in genes)) {
                name = (gene_name != "" ? gene_name : gene_id)
                # Track which requested identifiers were seen so we can report
                # any genes from the input list that are absent from the GTF.
                if (gene_id != "")   matched[gene_id] = 1
                if (gene_name != "") matched[gene_name] = 1

                print chrom, start, end, name
            }
        }

        END {
            # Write out only the requested genes that were actually matched in
            # the annotation so the shell can derive the missing subset.
            for (gene in genes) {
                if (gene in matched) {
                    print gene > "matched_genes.txt"
                }
            }
        }
    ' "${gtf}" | \
        sort -k1,1 -k2,2n -k3,3n -u > "${gtf.baseName}.genes.bed"
        # Sort and deduplicate the BED records to keep downstream interval
        # processing stable and deterministic.

    # Compare the requested genes with the matched genes and keep only the
    # identifiers that were never observed in the annotation.
    sort -u matched_genes.txt -o matched_genes.txt 2>/dev/null || touch matched_genes.txt
    comm -23 requested_genes.txt matched_genes.txt > missing_genes.txt

    # Fail early if none of the requested genes were found in the annotation.
    if [[ ! -s "${gtf.baseName}.genes.bed" ]]; then
        echo "ERROR: No matching genes from ${gene_list} were found in ${gtf}" >&2
        exit 1
    fi

    # For partial matches (subset of input genes are found in annotation), keep the pipeline running but point the user to the
    # file containing the genes in the input list that were absent from the annotation.
    if [[ ! -s missing_genes.txt ]]; then
        rm -f missing_genes.txt
    else
        echo "WARNING: Some genes in ${gene_list} were not found in ${gtf}. See missing_genes.txt for details." >&2
    fi
    """
}

// Compute panel on-target metrics from the merged gene-count matrix and panel gene list.
// The custom Python script writes MultiQC-compatible tabular output.
process CALC_PANEL_ON_TARGET {
    label "single"
    publishDir "${params.outdir}/panel_metrics", mode: 'copy'

    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/pandas:2.2.1'
        : 'quay.io/biocontainers/pandas:2.2.1'}"

    input:
    path gene_counts_matrix
    path panel_gene_list

    output:
    path "panel_metrics_on_target_mqc.tsv", emit: mqc_tsv

    script:
    """
    calc_panel_on_target.py -c ${gene_counts_matrix} -t ${panel_gene_list} -o panel_metrics_on_target_mqc.tsv
    """
}

// Calculate per-sample read count and N50 metrics for select genes.
process CALC_PANEL_SELECT_GENES_STATS {
    publishDir "${params.outdir}/panel_metrics", mode: 'copy'

    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.22.1--py39hdd5828d_3'
        : 'quay.io/biocontainers/pysam:0.22.1--py39hdd5828d_3'}"

    input:
    path bams
    path bais
    path bed

    output:
    path "panel_metrics_gene_stats_mqc.tsv", emit: mqc_tsv

    script:
    """
    calc_panel_select_genes_stats.py -b ${bams} -bed ${bed} -o panel_metrics_gene_stats_mqc.tsv
    """
}
