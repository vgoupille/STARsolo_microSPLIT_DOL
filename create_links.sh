#!/bin/bash

# This script creates the necessary directory structure and symbolic links
# for the STARsolo analysis script to find the input files.

# Base directories from the STARsolo script
DATA_DIR="raw_data"
FASTQ_DIR="${DATA_DIR}/fastq"
BARCODE_DIR="${DATA_DIR}/barcodes"
GENOME_REF_DIR="${DATA_DIR}/genome_ref"

# Source data locations
SOURCE_FASTQ="1_data/Fastq"
SOURCE_BARCODES="1_data/barcodes"
SOURCE_GENOME="1_data/genome_ref"

# Create directory structure
echo "Creating directory structure..."
mkdir -p $FASTQ_DIR
mkdir -p $BARCODE_DIR
mkdir -p $GENOME_REF_DIR

# Create symbolic links for FASTQ files
echo "Creating symbolic links for FASTQ files..."
ln -sf $(realpath $SOURCE_FASTQ/microSPLIT-600cells_S1_L001_R1_001.fastq.gz) $FASTQ_DIR/microSPLIT-600cells_S1_L001_R1_001.fastq.gz
ln -sf $(realpath $SOURCE_FASTQ/microSPLIT-600cells_S1_L001_R2_001.fastq.gz) $FASTQ_DIR/microSPLIT-600cells_S1_L001_R2_001.fastq.gz

# # Create symbolic links for barcode files
# echo "Creating symbolic links for barcode files..."
# ln -sf $(realpath $SOURCE_BARCODES/barcode_round1.txt) $BARCODE_DIR/barcode_round1.txt
# ln -sf $(realpath $SOURCE_BARCODES/barcode_round2.txt) $BARCODE_DIR/barcode_round2.txt
# ln -sf $(realpath $SOURCE_BARCODES/barcode_round3.txt) $BARCODE_DIR/barcode_round3.txt

# # Create symbolic links for genome reference files
# echo "Creating symbolic links for genome reference files..."
# ln -sf $(realpath $SOURCE_GENOME/GCA_030064105.1_ASM3006410v1_genomic.fna) $GENOME_REF_DIR/GCA_030064105.1_ASM3006410v1_genomic.fna
# ln -sf $(realpath $SOURCE_GENOME/genomic_annotation_R401.gff) $GENOME_REF_DIR/genomic_annotation_R401.gff

echo "Symbolic links created successfully."
echo "You can now run the STARsolo script which will find the input files at the expected locations." 