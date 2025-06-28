#!/usr/bin/env bash

# get test data from ONT's workflow: https://github.com/epi2me-labs/wf-transcriptomes

set -e

if [ -d "data/test/wf-transcriptomes-demo" ]
then
    echo "Directory already exists"
    exit 1
fi

# download data and uncompress
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-transcriptomes/wf-transcriptomes-demo.tar.gz
tar -xzvf wf-transcriptomes-demo.tar.gz

# remove archive
rm wf-transcriptomes-demo.tar.gz

# create directory for test data
mkdir -p data/test

# move the downloaded data into a temporary folder in the test data directory
# we download to temp dir so we can only the few files that are needed
# and get rid of the rest
mv wf-transcriptomes-demo/ data/test/wf-transcriptomes-demo-temp

# this directory will store files that are required
mkdir -p data/test/wf-transcriptomes-demo

# copy barcoded data folder
mv data/test/wf-transcriptomes-demo-temp/differential_expression_fastq/ data/test/wf-transcriptomes-demo/input_data

# copy genome files
mv data/test/wf-transcriptomes-demo-temp/hg38_chr20.fa data/test/wf-transcriptomes-demo/hg38_chr20.fa

mv data/test/wf-transcriptomes-demo-temp/gencode.v22.annotation.chr20.gtf data/test/wf-transcriptomes-demo/gencode.v22.annotation.chr20.gtf

# now that we have files we need, get rid of downloaded data
rm -r data/test/wf-transcriptomes-demo-temp
