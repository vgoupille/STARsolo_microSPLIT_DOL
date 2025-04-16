#!/bin/bash

# === Global Configuration for STARsolo Pipeline ===

# Here you can change the name of the analysis, the location of the folder for the Analysis, the number of threads, the memory, the maximum runtime, the location of the conda environment, the source data locations, the source file names, the target file names...

# Go to the folder where you want to run the analysis:
# For example: 
    # cd 1_data/test/starsolo_script_DOL_microSPLIT/
    # Then you can run the analysis with the following command:
    # sbatch run_pipeline.sh


# Email Configuration
EMAIL="valentin.goupille@univ-rennes.fr"

# Define the base path for all directories
BASE_PATH="/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/test/STARsolo_microSPLIT_DOL"  # <-- Change this to your actual base directory

# IMPORTANT: Make sure this path exists before running the pipeline
# You can check it using: [ -d "$BASE_PATH" ] && echo "Path exists" || echo "Path DOES NOT exist"
# If it doesn't exist, the pipeline will fail when trying to create the conda environment

# Analysis Name and Base Directory
ANALYSIS_NAME="STARsolo_microSPLIT" # For the name of the analysis: In the slurm output file, it will be the name of the job
BASE_DIR="Analysis_STARsolo_microSPLIT" # For the location of the folder for the Analysis: You can change it

# Resource Configuration
THREADS=16 # For the number of threads (cpus-per-task): You can change it # VERY IMPORTANT: YOU NEED TO HAVE ENOUGH THREADS TO RUN THE ANALYSIS: 16 MAY NOT BE ENOUGH FOR THE ANALYSIS =======> test with 32 or 64 is probably better

MEMORY="64G" # For the memory: You can change it # VERY IMPORTANT: YOU NEED TO HAVE ENOUGH MEMORY TO RUN THE ANALYSIS: 16G MAY NOT BE ENOUGH FOR THE ANALYSIS =======> test with 16G 32G or 64G IS PROBABLY BETTER

MAX_RUNTIME="10:00:00" # For the maximum runtime: You can change it

# Conda Configuration
CONDA_INIT_SCRIPT="/local/env/envconda.sh" # Path to the conda initialization script - Change this if your cluster uses a different location
CONDA_ENV_PATH="${BASE_PATH}/${BASE_DIR}/env_STARsolo" # For the location of the conda environment

# Source Data Locations - Absolute paths to prevent errors
SOURCE_FASTQ="/home/genouest/cnrs_umr6553/vgoupille/DOL_scRNAseq/1_data/Fastq"  # <-- Change this to your actual source data location


# for barcodes, genome_ref and genome_annotation here we use the file avaible in the repository (raw_data folder) :::::
SOURCE_BARCODES="${BASE_PATH}/raw_data/barcodes" 
SOURCE_GENOME_REF="${BASE_PATH}/raw_data/genome_ref" 
SOURCE_GENOME_ANNOTATION="${BASE_PATH}/raw_data/genome_annotation"





#---------------------------------------------------------------------------

# Source File Names - Update these to match your actual source files
SOURCE_FASTQ_R1="microSPLIT-600cells_S1_L001_R1_001.fastq.gz" # <-- Change this to your actual source file name
SOURCE_FASTQ_R2="microSPLIT-600cells_S1_L001_R2_001.fastq.gz" # <-- Change this to your actual source file name







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