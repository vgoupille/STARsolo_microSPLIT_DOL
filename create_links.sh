#!/bin/bash
#SBATCH --job-name=create_links_DOL_microsplit
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valentin.goupille@univ-rennes.fr

# This script creates the necessary directory structure and symbolic links
# for the STARsolo analysis script to find the input files.

# Base directories from the STARsolo script
BASE_DIR="Analysis_STARsolo_microsplit"
DATA_DIR="${BASE_DIR}/raw_data"
FASTQ_DIR="${DATA_DIR}/fastq"
BARCODE_DIR="${DATA_DIR}/barcodes"
GENOME_REF_DIR="${DATA_DIR}/genome_ref"
GENOME_ANNOTATION_DIR="${DATA_DIR}/genome_annotation"

# Source data locations - Update these to match your actual source directories
SOURCE_FASTQ="raw_data/fastq"
SOURCE_BARCODES="raw_data/barcodes"
SOURCE_GENOME_REF="raw_data/genome_ref"
SOURCE_GENOME_ANNOTATION="raw_data/genome_annotation"

# Create directory structure
echo "Creating directory structure..."
mkdir -p $FASTQ_DIR
mkdir -p $BARCODE_DIR
mkdir -p $GENOME_REF_DIR
mkdir -p $GENOME_ANNOTATION_DIR

# Create symbolic links for FASTQ files
echo "Creating symbolic links for FASTQ files..."
# Update these with your actual FASTQ filenames
ln -sf $(realpath $SOURCE_FASTQ/microSPLIT-600cells_S1_L001_R1_001.fastq.gz) $FASTQ_DIR/microSPLIT-600cells_S1_L001_R1_001.fastq.gz
ln -sf $(realpath $SOURCE_FASTQ/microSPLIT-600cells_S1_L001_R2_001.fastq.gz) $FASTQ_DIR/microSPLIT-600cells_S1_L001_R2_001.fastq.gz

# Create symbolic links for barcode files
echo "Creating symbolic links for barcode files..."
ln -sf $(realpath $SOURCE_BARCODES/barcode_round1.txt) $BARCODE_DIR/barcode_round1.txt
ln -sf $(realpath $SOURCE_BARCODES/barcode_round2.txt) $BARCODE_DIR/barcode_round2.txt
ln -sf $(realpath $SOURCE_BARCODES/barcode_round3.txt) $BARCODE_DIR/barcode_round3.txt

# Create symbolic links for genome reference files
echo "Creating symbolic links for genome reference files..."
ln -sf $(realpath $SOURCE_GENOME_REF/GCA_030064105.1_ASM3006410v1_genomic.fna) $GENOME_REF_DIR/GCA_030064105.1_ASM3006410v1_genomic.fna

# Create symbolic links for genome annotation files
echo "Creating symbolic links for genome annotation files..."
ln -sf $(realpath $SOURCE_GENOME_ANNOTATION/GCA_030064105.1_ASM3006410v1_genomic.gff) $GENOME_ANNOTATION_DIR/GCA_030064105.1_ASM3006410v1_genomic.gff

echo "Symbolic links created successfully."
echo "You can now run the STARsolo script which will find the input files at the expected locations." 