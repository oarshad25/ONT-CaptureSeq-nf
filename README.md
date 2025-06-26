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

The required parameters are as follows

| Parameter     | Type   | Description                                                                                                        | Default |
| ------------- | ------ | ------------------------------------------------------------------------------------------------------------------ | ------- |
| `inputdir`    | string | Path to parent read directory containing demultiplexed barcode subdirectories with fastq reads for each barcode OR |
| `samplesheet` | string | Path to sample manifest 'csv' file                                                                                 |
| `genome`      | string | Path to reference genome FASTA to use for alignment.                                                               |
| `annotation`  | string | Reference annotation as GTF.                                                                                       |

### Additional parameters

| Parameter              | Type    | Description                                                                                                                                                                               | Default   |
| ---------------------- | ------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------- |
| `outdir`               | string  | Output directory for results                                                                                                                                                              | `results` |
| `min_reads_per_sample` | integer | Threshold for minimum number of reads per sample. Samples with number of reads below this threshold are filtered out. This is mainly to get rid of unassigned barcodes in input directory | 1000      |

