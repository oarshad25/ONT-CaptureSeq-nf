#!/usr/bin/env nextflow

/*
 * Drop FASTQ files based on read count
    * This subworkflow takes input FASTQ files and drops those that have a read count below a specified threshold.
    * The samples that are dropped if there are any are logged in a CSV file with columns for sample ID, filename, and read count.
*/

// use seqkit stats to add text file with read count to input fastq channel
process COUNT_READS {
    tag "${meta.id}"
    label "very_short"

    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    input:
    // input fastq
    tuple val(meta), path(fastq)

    output:
    // output the original input channel with an additional path to a text file containing the read count for the fastq
    tuple val(meta), path(fastq), path("count.txt")

    script:
    """
    # use seqkit stats to get stats on fastq, skip header (NR==2) and print 4th column (num_seqs) to file
    seqkit stats \\
        --tabular \\
        --threads ${task.cpus} \\
        ${fastq} | \\
    awk 'NR==2 {print \$4}' > count.txt
    """
}

workflow DROP_FASTQS_ON_READ_COUNT {
    take:
    reads_ch // tuple val(meta), path(fastq); queue channel with tuple of metamap and fastq
    min_reads // integer. Any samples (fastq's) in reads_ch with read count below this threshold will be dropped
    step // string, used to separate different invocations of this workflow such as on input fastqs vs filtered fastqs. The location of log file is determined by this

    main:
    // Count number of reads in each input FASTQ file using seqkit stats and add this as a text file to the input channel
    COUNT_READS(reads_ch)

    // Branch samples into a "passed" and "failed" channel based on whether they meet the read count threshold
    COUNT_READS.out
        .branch { _meta, _fastq, count_file ->
            // Read the read count from the count file
            // .text is syntactic sugar for .getText()
            def count = count_file.text.trim().toInteger()
            passed: count >= min_reads
            failed: count < min_reads
        }
        .set { read_counts_ch_branched }

    // Reshape the passed reads counts channel to drop the count file to match input channel caardinality
    read_counts_ch_branched.passed
        .map { meta, fastq, _count_file ->
            tuple(meta, fastq)
        }
        .set { passed_reads_ch }

    // Write failed samples to a CSV file if there are any, otherwise skip

    // collect all failed items into a list so that we can check if there any.
    // Cannot use if, channels are objects and don't evaluate to true/false, also cannot use .size() on a streaming channel
    read_counts_ch_branched.failed
        .toList()
        .filter { failed_list ->
            // Only continue down this stream if the list has at least 1 item
            failed_list.size() > 0
        }
        .map { failed_list ->
            // Start building the CSV content with the header
            def csv_content = "Sample_ID,Filename,Read_Count\n"

            // Loop through the failed items and append their data
            failed_list.each { meta, fastq, count_file ->
                def count = count_file.text.trim()
                csv_content += "${meta.id},${fastq.name},${count}\n"
            }

            return csv_content
        }
        .collectFile(
            name: 'dropped_samples.csv',
            storeDir: "${params.outdir}/qc/read_count_filtering/${step}",
            newLine: false,
        )
        .view { "📝 Wrote dropped sample log (step: ${step}) to: ${it}" }

    emit:
    passed_reads_ch // tuple val(meta), path(fastq) of samples that passed read count threshold
}
