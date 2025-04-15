#!/bin/bash

# === Global Configuration for STARsolo Pipeline ===

#Here you can change the name of the analysis, the location of the folder for the Analysis, the number of threads, the memory, the maximum runtime, the location of the conda environment, the source data locations, the source file names, the target file names...

#Go to the folder where you want to run the analysis :
# for example : 
    #cd 1_data/test/starsolo_script_DOL_microSPLIT/
    #then you can run the analysis with the following command :
    #sbatch run_pipeline.sh


# Email Configuration
EMAIL="valentin.goupille@univ-rennes.fr"

# Analysis Name and Base Directory
ANALYSIS_NAME="STARsolo_DOL_microSPLIT" #for the name of the analysis : In the slurm output file, it will be the name of the job
BASE_DIR="Analysis_STARsolo_microSPLIT" #for the location of the folder for the Analysis : You can change it

# Resource Configuration
THREADS=16 #for the number of threads (cpus-per-task): You can change it #VERY IMPORTANT : YOU NEED TO HAVE ENOUGH THREADS TO RUN THE ANALYSIS : 16 IS MABY NOT ENOUGH FOR THE ANALYSIS =======> test with 32 or 64 is probably better

MEMORY="16G" #for the memory : You can change it #VERY IMPORTANT : YOU NEED TO HAVE ENOUGH MEMORY TO RUN THE ANALYSIS : 16G IS MABY NOT ENOUGH FOR THE ANALYSIS =======> test with 32G or 64G IS PROBABLY BETTER

MAX_RUNTIME="10:00:00" #for the maximum runtime : You can change it

# Environment Paths
CONDA_ENV_PATH="${BASE_DIR}/env_STARsolo" #for the location of the conda environment : You can change it (not necessary in Analysis_STARsolo_microSPLIT)

# Source Data Locations - Absolute paths to prevent errors
# Define your base path here (change this to your actual base path)
BASE_PATH="/path/to/your/base/directory"  # <-- Change this to your actual base directory

# Access to files in the current directory
SOURCE_FASTQ="/path/to/your/fastq/directory" 
SOURCE_BARCODES="${BASE_PATH}/barcodes/directory" 

# Access to files in the parent directory (one level up)
SOURCE_GENOME_REF="${BASE_PATH}/../shared_references/genome_ref/directory" 
# Example with direct absolute path without using variables
SOURCE_GENOME_ANNOTATION="/absolute/path/to/genome_annotation/directory"

# Example of accessing a file two levels up
# ANOTHER_REFERENCE="${BASE_PATH}/../../global_references/other_file"

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