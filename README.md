# ONT-CaptureSeq-nf

**ONT-CaptureSeq-nf** is a bioinformatics pipeline for the analysis of
multiplexed long-read nanopore targeted RNA capture sequencing developed in [Nextflow](https://www.nextflow.io/).

## Getting ONT-CaptureSeq-nf
```bash
$ git clone git@github.com:oarshad25/ONT-CaptureSeq-nf.git
$ cd ONT-CaptureSeq-nf
```
## Usage

```bash
$ nextflow run main.nf --help
```

Example usage:

```bash
$ nextflow run main.nf \
  --inputdir "path/to/inputdir" \
  --genome "path/to/genome.fa" \
  --annotation "path/to/annotation.gtf"
```

Can also use params.json to specify the parameters

```bash
$ nextflow run main.nf -params-file params.json
```

### Input

The workflow accepts either an input directory (specified via `--inputdir`)
or a sample manifest (specified via `--samplesheet`).

(a) *input directory:* Path to demultiplexed barcode directory, the directory
containing one level of sub-directories which in turn contain FASTQ files

```
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

(b) *sample manifest:* Path to sample sheet in csv format. The sample manifest should be a comma seperated values (.csv)
file and include at least two columns named `id` and `fastqdir`

- the `id` column is the sample id for each sample/barcode
- the `fastqdir` column is the path to the directory containing the fastq files for the barcode.


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
| `min_qual`             | float   | Average read quality threshold. Reads below this threshold are filtered out                                                                                                               | 7       |

#### Restrander (read orientation)

| Parameter           | Type    | Description                                | Default                                |
| ------------------- | ------- | ------------------------------------------ | -------------------------------------- |
| `run_restrander`    | boolean | Reorient reads with restrander             | false                                  |
| `restrander_config` | string  | Path to Restrander configuration json file | `assets/restrander_config/PCB109.json` |

#### Alignment

##### MiniMap2

Command-line options for Minimap2 alignment. See [minimap2 options](https://lh3.github.io/minimap2/minimap2.html)

| Parameter                      | Type      | Description                                                                           | Default                     |
| ------------------------------ | --------- | ------------------------------------------------------------------------------------- | --------------------------- |
| `skip_save_minimap2_index`     | boolean   | Whether to skip prebuilding and saving Minimap2 index.                                | false                       |
| `minimap2_indexing_extra_opts` | string    | Any extra options to be provided to Minimap2 indexing                                 |                             |
| `minimap2_junc_bed`            | string    | Optional path to Minimap junction bed file (used by `--junc-bed`)                     |                             |
| `minimap2_x`                   | string    | Preset, Minimap2 `-x` option e.g. `"map-ont"`                                         | `"splice"`                  |
| `minimap2_k`                   | integer   | kmer size, Minimap2 `-k` option e.g. `14`                                             | Minimap2 default (`15`)     |
| `minimap2_u`                   | character | Minimap2 `-u` option e.g. `"b"`                                                       | Minimap2 default (`"n"`)    |
| `minimap2_G`                   | string    | Maximum intron length, Minimap2 `-G` option e.g. `"500K"`                             | Minimap2 default (`"200k"`) |
| `minimap2_I`                   | string    | Batchsize for indexing, Minimap2 `-I` option e.g. `"8G"`                              | Minimap2 default            |
| `minimap2_cs`                  | string    | Output cs tag, Minimap2 `--cs` option e.g. `"long"`                                   | Minimap2 default (none)     |
| `minimap2_extra_opts`          | string    | Any extra options to be provided to Minimap2 e.g. `"--splice-flank=no"` for SIRV data |                             |

##### RSeQC

| Parameter    | Type    | Description                                                    | Default |
| ------------ | ------- | -------------------------------------------------------------- | ------- |
| `skip_rseqc` | boolean | skip read distribution calculation of aligned reads with RSeQC | false   |

##### Additional

| Parameter           | Type    | Description                                                                                                  | Default |
| ------------------- | ------- | ------------------------------------------------------------------------------------------------------------ | ------- |
| `filter_bam_mapped` | boolean | Whether to filter aligned bam to mapped reads only i.e. whether to filter out unmapped reads from alignments | true    |

#### Isoform Discovery and Quantification

##### IsoQuant

| Parameter                   | Type    | Description                                                                                         | Default |
| ----------------------------| ------- | ----------------------------------------------------------------------------------------------------| ------- |
| `isoquant_complete_genedb`  | boolean | Whether to set option --complete_genedb. Set to true for official annotations such as from GENCODE  | true    |