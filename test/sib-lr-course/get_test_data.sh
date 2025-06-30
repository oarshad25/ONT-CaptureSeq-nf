#!/usr/bin/env bash

# get test data from SIB long read analysis course version 2025.2 project1:
# https://sib-swiss.github.io/NGS-longreads-training/2025.2/course_material/group_work/project1/

set -e

shopt -s extglob

TEST_DATA_DIR="data/test/sib-lr-course/"

if [ -d $TEST_DATA_DIR ]
then
    echo "Test data directory $TEST_DATA_DIR already exists"
    exit 1
fi

# download data
wget https://ngs-longreads-training.s3.eu-central-1.amazonaws.com/project1.tar.gz
tar -xvf project1.tar.gz
rm project1.tar.gz

# create a directory for the fastq reads
READ_DIR="${TEST_DATA_DIR}reads"
mkdir -p $READ_DIR

# move reads for each sample into it's own subdirectory within 'READ_DIR'
for fastq in project1/reads/*.fastq.gz
do
    sample_name=$(basename "$fastq" .fastq.gz)
    # get rid of '_'
    sample_name=${sample_name//_/}

    mkdir -p "${READ_DIR}/${sample_name}"
    cp $fastq "${READ_DIR}/${sample_name}/${sample_name}.fastq.gz"
done

# copy references
mkdir "${TEST_DATA_DIR}references"
cp -r project1/references "${TEST_DATA_DIR}references"

rm -r project1

