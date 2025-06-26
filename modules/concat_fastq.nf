#!/usr/bin/env nextflow

/*
================================================================================
    CONCAT_FASTQ module
================================================================================

Concatenate/merge list of fastq files

Take a channel of the shape `[meta, fastqlist]` and merge the files in fastqlist.
Output a channel with merged fastq file added to the metamap

================================================================================
*/

/*
name: concat_fastq
description: Concatenates fastq files
tools:
  - cat:
      description: |
        The cat utility reads files sequentially, writing them to the standard output.
        documentation: https://www.gnu.org/software/coreutils/manual/html_node/cat-invocation.html
      licence: ["GPL-3.0-or-later"]
      identifier: ""
input:
  - - meta:
        type: map
        description: |
          Groovy Map containing sample information
          e.g. [ id:'barcode01' ]
    - fastq_files:
        type: list of files
        description: |
          List of input FastQ files (optionally compressed) to be concatenated.
          Files must end in extension .fastq (or fastq.gz)

output:
  - merged_reads:
      - meta:
          type: map
          description: |
            Groovy Map containing sample information
            e.g. [ id:'test' ]
      - "*.merged.fastq(.gz)":
          type: file
          description: Merged fastq file
*/

process CONCAT_FASTQ {
    tag { meta.id }

    publishDir "${params.outdir}/merged_reads", mode: 'copy'

    input:
    // Reads to merge
    // [meta, [ read1, read2 ] ]
    tuple val(meta), path(fastq_files)

    output:
    // merged fastq file
    tuple val(meta), path("${meta.id}.merged.${extn}"), emit: merged_reads

    script:
    if (fastq_files.every { it.name.endsWith('.fastq.gz') }) {
        extn = 'fastq.gz'
    }
    else if (fastq_files.every { it.name.endsWith('.fastq') }) {
        extn = 'fastq'
    }
    else {
        error("Concatentation of mixed filetypes is unsupported")
    }

    """
    cat ${fastq_files} > "${meta.id}.merged.${extn}"
    """
}
