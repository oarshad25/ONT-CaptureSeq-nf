# ONT-CaptureSeq-nf

**ONT-CaptureSeq-nf** is a bioinformatics pipeline for the analysis of
multiplexed long-read nanopore targeted RNA capture sequencing.

The pipeline is developed in [Nextflow](https://www.nextflow.io/) with software
dependencies managed using [Docker](https://www.docker.com) or [Apptainer](https://apptainer.org).

## Pipeline Summary

1. Prepare reference files
    * Unzip reference fasta and GTF if compressed
    * Convert GTF to BED using UCSC utilities [ucsc-gtftogenepred](https://open.bioqueue.org/home/knowledge/showKnowledge/sig/ucsc-gtftogenepred)
      and [ucsc-genepredtobed](https://open.bioqueue.org/home/knowledge/showKnowledge/sig/ucsc-genepredtobed).
2. Merge fastq files per sample and decompress  ([zcat](https://linux.die.net/man/1/zcat)).
3. Drop any samples with small number of reads by filtering samples with read count below a threshold.
4. Raw read QC ([NanoPlot](https://github.com/wdecoster/NanoPlot), [NanoComp](https://github.com/wdecoster/nanocomp) [MultiQC](https://multiqc.info/docs/)).
5. Filter reads on length and quality ([seqkit seq](https://bioinf.shenwei.me/seqkit/usage/#seq); _optional_).
    * Summary statistics of filtered reads.
6. Preprocessing of cDNA reads with one of the following (_optional_):
    * Identify, orient and trim full-length reads ([pychopper](https://github.com/epi2me-labs/pychopper))
    * Reorient reads ([restrander](https://github.com/mritchielab/restrander)).
   Generate summary statistics on QC'ed reads with ([seqkit stats](https://bioinf.shenwei.me/seqkit/usage/#stats))
7. Align reads to genome ([minimap2](https://github.com/lh3/minimap2))
    * Optionally prebuild minimap2 genome index.
    * Optionally use reference annotation to prioritise on known splice junctions.
8. Sort, index alignments and generate alignment statistics ([samtools](https://www.htslib.org/doc/))
9. Additional QC of aligned reads
    * [Cramino](https://github.com/wdecoster/cramino).
    * Optionally [RSeQC](https://rseqc.sourceforge.net).
10. Summarise QC stats for processed fastq's and mapped reads ([MultiQC](https://multiqc.info/docs/)).
11. Filter out unmapped, secondary and supplementary reads ([samtools](https://www.htslib.org/doc/)).
12. Optionally, if a genelist of interest is provided, generate a subset of aligned reads to genes in the list.
    * QC metrics of subset of aligned reads to assess quality of reads mapping to genes of interest (skipped for pychopper processed reads).
13. Generate gene level count matrix with [featureCounts](https://subread.sourceforge.net/featureCounts.html).
14. Transcript reconstruction and quantification ([IsoQuant](https://ablab.github.io/IsoQuant/) or [FLAIR](https://flair.readthedocs.io/en/latest/index.html)).

## Getting ONT-CaptureSeq-nf

```bash
git clone git@github.com:oarshad25/ONT-CaptureSeq-nf.git
cd ONT-CaptureSeq-nf
```

## Usage

```bash
nextflow run main.nf --help
```

Example usage:

```bash
nextflow run main.nf \
  -profile <docker/apptainer> \
  --samplesheet "path/to/samplesheet" \
  --genome "path/to/genome.fa" \
  --annotation "path/to/annotation.gtf"
```

In addition to specifying parameters at the command line, can also specify the
parameters through either editing the parameter configuration file [./conf/params.config](conf/params.config)
or via Nextflow `-params-file` option:

```bash
nextflow run main.nf -profile <docker/apptainer> -params-file params.json
```

Can customise execution to the specific compute platform on which the pipeline is executed, by editing the `nextflow.config` file.
This pipeline was developed on University of Oxford's [BMRC](https://www.medsci.ox.ac.uk/for-staff/resources/bmrc) cluster:

```bash
$ nextflow run main.nf \
  -profile bmrc \
  --samplesheet "path/to/samplesheet" \
  --genome "path/to/genome.fa" \
  --annotation "path/to/annotation.gtf"
```

### Input

The workflow accepts either an input directory (specified via `--inputdir`)
or a sample manifest (specified via `--samplesheet`).

(a) _input directory:_ Path to demultiplexed barcode directory, the directory
containing one level of sub-directories which in turn contain FASTQ files

``` bash
─── inputdir
    ├── barcode01
    │   ├── reads0.fastq
    │   └── reads1.fastq
    ├── barcode02
    │   ├── reads0.fastq
    │   ├── reads1.fastq
    │   └── reads2.fastq
    └── barcode03
        └── reads0.fastq
```

OR

(b) _sample manifest:_ Path to sample sheet in csv format. The sample manifest should be a comma seperated values (.csv)
file and include at least two columns named `id` and `fastqdir`

* the `id` column is the sample id for each sample/barcode
* the `fastqdir` column is the path to the directory containing the fastq files for the barcode.

### Pipeline parameters

The required parameters are as follows:

| Parameter     | Type   | Description                                                                                                        | Default |
| ------------- | ------ | ------------------------------------------------------------------------------------------------------------------ | ------- |
| `inputdir`    | string | Path to parent read directory containing demultiplexed barcode subdirectories with fastq reads for each barcode OR |
| `samplesheet` | string | Path to sample manifest 'csv' file                                                                                 |
| `genome`      | string | Path to reference genome FASTA to use for alignment.                                                               |
| `annotation`  | string | Reference annotation as GTF.                                                                                       |

### Optional parameters

| Parameter | Type   | Description                  | Default |
| --------- | ------ | ---------------------------- | ------- |
| `outdir`  | string | Output directory for results | results |

### Additional Options

Options for configuring steps/tools in the workflow

#### QC and filtering

| Parameter              | Type    | Description                                                                                                                                                                               | Default |
| ---------------------- | ------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------- |
| `min_reads_per_sample` | integer | Threshold for minimum number of reads per sample. Samples with number of reads below this threshold are filtered out. This is mainly to get rid of unassigned barcodes in input directory | 1000    |
| `is_fastq_rich`        | boolean | Used in NanoPlot. Whether input is a rich fastq with additional information regarding ONT run (concerning channel and time). Used to set input data source argument for nanoplot          | false   |
| `skip_read_filtering`  | boolean | Skip filtering of reads based on length and average quality                                                                                                                               | false   |
| `min_length`           | integer | Minimum read length threshold. Reads below this length are filtered out                                                                                                                   | 100     |
| `max_length`           | integer | Maximum read length threshold. Reads longer than this length are filtered out                                                                                                             | 25000   |
| `min_qual`             | float   | Average read quality threshold. Reads below this threshold are filtered out                                                                                                               | 7       |

##### cDNA reads preprocessing

Options for further preprocessing of cDNA reads.

| Parameter           | Type    | Description                                | Default                                |
| ------------------- | ------- | ------------------------------------------ | -------------------------------------- |
| `run_cdna_qc` | boolean       | Whether to preprocess reads using a cDNA readQC tool | true |
| `cdna_qc_method` | boolean       | Which tool to use for cDNA read QC/preprocessing. Options are `"pychopper"` or `"restrander"`. | `pychopper` |

###### Pychopper (full length reads)

| Parameter           | Type    | Description                                | Default                                |
| ------------------- | ------- | ------------------------------------------ | -------------------------------------- |
| `pychopper_cdna_kit` | string  | kit name, pychopper -k option e.g. "PCS109" |  |
| `pychopper_primer_fasta` | string | Path to primer sequences fasta, pychopper -b option |  |
| `pychopper_backend` | string | Backend to be used by pychopper for identifying primers in input reads, pychopper -m option (phmm or edlib) | `"edlib"`|
| `pychopper_extra_opts` | string | Any extra command line options to pass to pychopper e.g. `"-p"` to keep primers | |

If neither `pychopper_cdna_kit` nor `pychopper_primer_fasta` is provided, pychopper default (kit PCS109) is used.

###### Restrander (read orientation)

| Parameter           | Type    | Description                                | Default                                |
| ------------------- | ------- | ------------------------------------------ | -------------------------------------- |
| `restrander_config` | string  | Path to Restrander configuration json file | `assets/restrander_config/PCB109.json` |

#### Alignment

##### General Options

| Parameter                       | Type    | Description                                                                                                                                                                                   | Default |
| ------------------------------- | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------- |
| `alignment_use_annotation`      | boolean | Minimap2 can optionally take annotated genes as input. This parameter selects whether to use reference annotation in Minimap2 alignment as input to prioritise on annotated splice junctions. | false   |
| `filter_alignments`             | boolean | Whether to filter alignment bams to remove secondary, supplementary and unmapped reads                                                                                                        | true    |
| `genelist_to_subset_alignments` | string  | Path to optional file containing list of gene id's, one per line. If this file is provided, a subset of the filtered aligned bams is generated that overlaps genes in this list               |         |

##### MiniMap2

Command-line options for Minimap2 alignment. See [minimap2 options](https://lh3.github.io/minimap2/minimap2.html)

| Parameter                      | Type      | Description                                                                                                                                                                               | Default                     |
| ------------------------------ | --------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------- |
| `skip_save_minimap2_index`     | boolean   | Whether to skip prebuilding and saving Minimap2 index.                                                                                                                                    | false                       |
| `minimap2_indexing_extra_opts` | string    | Any extra options to be provided to Minimap2 indexing                                                                                                                                     |                             |
| `minimap2_junc_bed`            | string    | Optional path to Minimap junction bed file (used by `--junc-bed`). Cannot be used in conjunction with `alignment_use_annotation` as that sets input annotation to be used by `--junc-bed` |                             |
| `minimap2_x`                   | string    | Preset, Minimap2 `-x` option e.g. `"map-ont"`                                                                                                                                             | `"splice"`                  |
| `minimap2_k`                   | integer   | kmer size, Minimap2 `-k` option e.g. `14`                                                                                                                                                 | Minimap2 default (`15`)     |
| `minimap2_u`                   | character | Minimap2 `-u` option e.g. `"b"`                                                                                                                                                           | Minimap2 default (`"n"`)    |
| `minimap2_G`                   | string    | Maximum intron length, Minimap2 `-G` option e.g. `"500K"`                                                                                                                                 | Minimap2 default (`"200k"`) |
| `minimap2_I`                   | string    | Batchsize for indexing, Minimap2 `-I` option e.g. `"8G"`                                                                                                                                  | Minimap2 default            |
| `minimap2_cs`                  | string    | Output cs tag, Minimap2 `--cs` option e.g. `"long"`                                                                                                                                       | Minimap2 default (none)     |
| `minimap2_junc_bonus`          | integer   | Score bonus for splice site in annotation. Effective if `alignment_use_annotation` is set or `minimap2_junc_bed` is provided                                                              | Minimap2 default (9)        |
| `minimap2_extra_opts`          | string    | Any extra options to be provided to Minimap2 e.g. `"--splice-flank=no"` for SIRV data                                                                                                     |                             |

##### RSeQC

| Parameter    | Type    | Description             | Default |
| ------------ | ------- | ----------------------- | ------- |
| `skip_rseqc` | boolean | skip read QC with RSeQC | false   |

#### Isoform Discovery and Quantification

| Parameter                  | Type    | Description                                                                       | Default      |
| -------------------------- | ------- | --------------------------------------------------------------------------------- | ------------ |
| `skip_isoform_discovery`   | boolean | Whether to skip running isoform discovery subworkflow                             | false        |
| `isoform_discovery_method` | string  | Which method to use for isoform discovery. Options are `"isoquant"` or `"flair"`. | `"isoquant"` |

##### IsoQuant

| Parameter                  | Type    | Description                                                                                                                    | Default |
| -------------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------ | ------- |
| `isoquant_complete_genedb` | boolean | Whether to set option --complete_genedb. Set to true for official annotations such as from GENCODE                             | true    |
| `isoquant_extra_opts`      | string  | Any additional command line options to pass to IsoQuant e.g. `"--sqanti_output --splice_correction_strategy conservative_ont"` |         |

##### Flair

| Parameter                    | Type   | Description                                                            | Default                                                                     |
| ---------------------------- | ------ | ---------------------------------------------------------------------- | --------------------------------------------------------------------------- |
| `flair_collapse_extra_opts`  | string | Any additional command line options to pass to flair collapse module   | `"--stringent --check_splice --generate_map --annotation_reliant generate"` |
| `flair_align_reads_manifest` | string | Required to run flair. Path to sample manifest for flair align module. |                                                                             |
