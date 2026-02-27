#!/usr/bin/env nextflow

/*
    * Workflow to drop fastq files that have a read count below a specified threshold.
    * This is used to filter out unassigned barcodes or samples that have very low read counts
    * The workflow also generates a CSV file listing the samples that were dropped due to low read count, along with their sample ID, filename and read count.
    * If there are no samples that fail the read count filter, the CSV file is not generated.
    * The workflow also generates seqkit stats files for the samples that pass the read count filter,
    * which are output to the main workflow to be used by the multiqc module.
    *
    * Inputs:
    * - reads_ch: channel with tuple of metamap and path to fastq file
    * - min_reads: integer. Any samples (fastq's) in reads_ch with read count below this threshold will be dropped
    * - step: string, parameter used by SEQKIT_STATS module to seperate mutliple invocations of seqkit stats
    *         e.g. on raw input fastqs vs filtered fastqs. This is used to give different names to the stats files
    *         and generate seperate tabs in the main multiqc report through the config
    *
    * Outputs:
    * - passed_reads_ch: channel with tuple of metamap and path to fastq file for samples that pass the read count filter
    * - seqkit_stats_ch: channel with tuple of metamap and path to seqkit stats tsv file for samples that pass the read count filter.
    *                    These stats files are used in multiqc report on raw reads that pass read count filter.
*/

include { SEQKIT_STATS } from '../modules/seqkit_stats.nf'

workflow DROP_FASTQS_ON_READ_COUNT {
    take:
    reads_ch // tuple val(meta), path(fastq); queue channel with tuple of metamap and fastq
    min_reads // integer. Any samples (fastq's) in reads_ch with read count below this threshold will be dropped
    step // string, used to separate different invocations of this workflow such as on raw input fastqs vs filtered fastqs. The location of log file is determined by this

    main:
    // 1) Generate stats on input reads

    // generate statistitcs on input reads using seqkit stats
    // these are used to filter out samples with low read counts before downstream processing
    // and also to generate stats for multiqc report on raw reads that pass read count filter
    SEQKIT_STATS(reads_ch, step)

    // 2) Add read count to seqkit stats output channel, combine with reads channel
    //    and branch into "passed" and "failed" channels based on whether samples meet the read count threshold

    // * From the stats file derive the read count by reading the text file with .text which is syntactic sugar for getText(),
    //   getting second row (skipping header) and pulling out the 4th column (num_seqs) which is the read count
    // * Join the read count with the channel containing the merged fastq files
    // * Filter out samples with read count less than `min_reads`
    // * Branch samples into a "passed" and "failed" channel based on whether they meet the read count threshold
    SEQKIT_STATS.out.stats
        .map { meta, tsv ->
            def count = tsv.text.split('\n')[1].split('\t')[3] as Integer
            [meta, tsv, count]
        }
        .join(reads_ch)
        .branch { _meta, _tsv, count, _fastq ->
            passed: count >= min_reads
            failed: count < min_reads
        }
        .set { branched_reads_and_stats_ch }

    // 3) For samples that pass read count filter, emit the fastq and stats file in separate channels
    //    of approriate shape and cardinality for downstream processing and multiqc report generation respectively.
    //    These channels form the output of this workflow.

    // MultiMap the filtered channel to get separate channels of requisite caardinality for the seqkit stats
    // files and the fastq files for the samples that pass the read count filter.
    // These are used in downstream steps and also to generate multiqc report on raw reads that pass read count filter.
    branched_reads_and_stats_ch.passed
        .multiMap { meta, tsv, _count, fastq ->
            stats: tuple(meta, tsv)
            reads: tuple(meta, fastq)
        }
        .set { passed_reads_and_stats_ch }

    passed_reads_ch = passed_reads_and_stats_ch.reads
    seqkit_stats_ch = passed_reads_and_stats_ch.stats

    // 4) For samples that fail read count filter, write out a CSV file with sample ID, filename and read count for each failed sample.
    //    If there are no failed samples, skip this step.

    // Write failed samples to a CSV file if there are any, otherwise skip
    // collect all failed items into a list so that we can check if there any.
    // Cannot use if, channels are objects and don't evaluate to true/false, also cannot use .size() on a streaming channel
    branched_reads_and_stats_ch.failed
        .toList()
        .filter { failed_list ->
            // Only continue down this stream if the list has at least 1 item
            failed_list.size() > 0
        }
        .map { failed_list ->
            // Start building the CSV content with the header
            def csv_content = "Sample_ID,Filename,Read_Count\n"

            // Loop through the failed items and append their data
            failed_list.each { meta, _tsv, count, fastq ->
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
    seqkit_stats_ch // tuple val(meta), path(stats_file) of samples that passed read count threshold
}
