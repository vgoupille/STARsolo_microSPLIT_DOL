#!/bin/bash

# === Global Configuration for STARsolo Pipeline ===

# Email Configuration
EMAIL="valentin.goupille@univ-rennes.fr"

# Analysis Name and Base Directory
ANALYSIS_NAME="STARsolo_DOL_microSPLIT"
BASE_DIR="Analysis_STARsolo_microsplit"

# Resource Configuration
THREADS=16
MEMORY="16G"
MAX_RUNTIME="4:00:00"

# Environment Paths
CONDA_ENV_PATH="/DOL_scRNAseq/5_environnements/env_STARsolo"

# Source Data Locations - Update these to match your actual source directories
SOURCE_FASTQ="1_data/Fastq" 
SOURCE_BARCODES="starsolo_script_DOL_microsplit/raw_data/barcodes" 
SOURCE_GENOME_REF="starsolo_script_DOL_microsplit/raw_data/genome_ref" 
SOURCE_GENOME_ANNOTATION="starsolo_script_DOL_microsplit/raw_data/genome_annotation"

# Source File Names - Update these to match your actual source files
SOURCE_FASTQ_R1="microSPLIT-600cells_S1_L001_R1_001.fastq.gz"
SOURCE_FASTQ_R2="microSPLIT-600cells_S1_L001_R2_001.fastq.gz"
SOURCE_GENOME_FASTA="GCA_030064105.1_ASM3006410v1_genomic.fna"
SOURCE_GENOME_GFF="GCA_030064105.1_ASM3006410v1_genomic.gff"

# Target (Symbolic Link) File Names
TARGET_FASTQ_R1="microSPLIT-R1.fastq.gz"
TARGET_FASTQ_R2="microSPLIT-R2.fastq.gz"
TARGET_GENOME_FASTA="genome_ref_PsR401.fna"
TARGET_GENOME_GFF="genome_annotation_PsR401.gff"

# Path Structure - Don't change these unless necessary
DATA_DIR="${BASE_DIR}/raw_data"
FASTQ_DIR="${DATA_DIR}/fastq"
BARCODE_DIR="${DATA_DIR}/barcodes"
GENOME_REF_DIR="${DATA_DIR}/genome_ref"
GENOME_ANNOTATION_DIR="${DATA_DIR}/genome_annotation"
GENERAL_OUTPUT_DIR="${BASE_DIR}/Output_DOL_microsplit_starsolo"
GENOME_DIR="${GENERAL_OUTPUT_DIR}/genome_index"
OUTPUT_STARSOLO_DIR="${GENERAL_OUTPUT_DIR}/starsolo_output" 