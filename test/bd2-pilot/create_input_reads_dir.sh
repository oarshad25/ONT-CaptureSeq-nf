#!/usr/bin/env bash

# this script copies a subset of barcodes for the BD2 pilot to create a test reads directory

# basecalled and demultiplexed reads are located in $FROM_DIR
# copy the reads for selected barcodes to the destination directory $TO_DIR

set -e

# https://stackoverflow.com/questions/19115156/show-commands-without-executing-them
debug=
#debug=echo

# set of barcodes to copy: array of two digit barcode numbers
barcodes=(59 61)
# echo "${barcodes[@]}"

# path to parent directory containing demultiplexed reads for BD2 pilot (this has subdirectories for each barcode)
FROM_DIR="/well/pharrison/users/piz649/projects/bd2_pilot/analysis/01_primary/demux_reads/results/output_reads"
# path to parent directory where reads for selected barcodes are to be copied to
TO_DIR="/well/pharrison/users/piz649/projects/ONT-CaptureSeq-nf/data/test/bd2-pilot/reads"

for barcode in "${barcodes[@]}"
do
   $debug cp -r "${FROM_DIR}/barcode${barcode}" "${TO_DIR}/barcode${barcode}"
done

