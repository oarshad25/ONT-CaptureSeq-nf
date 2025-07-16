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
  - zcat:
      description: |
        Concatenate (compressed) files
        documentation: https://linux.die.net/man/1/zcat
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
      - "*.fastq":
          type: file
          description: Merged fastq file
*/

process CONCAT_FASTQ {
  tag { meta.id }
  label "single"
  label "local_software"

  publishDir "${params.outdir}/merged_reads", mode: 'link'

  container "${workflow.containerEngine == 'apptainer'
    ? 'https://depot.galaxyproject.org/singularity/ubuntu:24.04'
    : 'quay.io/biocontainers/ubuntu:24.04'}"

  input:
  // Reads to merge
  // [meta, [ read1, read2 ] ]
  tuple val(meta), path(fastq_files)

  output:
  // merged fastq file
  tuple val(meta), path("${meta.id}.fastq"), emit: merged_reads

  script:
  """
  zcat -f ${fastq_files} > "${meta.id}.fastq"
  """
}
