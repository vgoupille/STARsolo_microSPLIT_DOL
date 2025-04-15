#!/bin/bash

# SLURM Job Configuration
#SBATCH --job-name=STARsolo_DOL_microSPLIT  # Job name
#SBATCH --mail-user=valentin.goupille@univ-rennes.fr
#SBATCH --mail-type=BEGIN,END,FAIL          # Email notifications
#SBATCH --ntasks=1                   # Number of tasks
#SBATCH --cpus-per-task=16           # CPUs per task
#SBATCH --mem=64G                    # Memory allocation
#SBATCH --time=24:00:00              # Maximum job runtime
#SBATCH --output=starsolo_%j.out     # Standard output log file
#SBATCH --error=starsolo_%j.err      # Standard error log file

# Conda Environment Setup
# Source the Conda initialization script
. /local/env/envconda.sh

# Activate the specific Conda environment for STARsolo
conda activate env_STARsolo  # Path to STARsolo environment

# Function to check command success
check_success() {
    if [ $? -ne 0 ]; then
        echo "ERROR: Command failed - $1"
        exit 1
    fi
}

# Global Configuration Variables
# Set number of threads for parallel processing
THREADS=16

# ===================================================================
# INPUT/OUTPUT PATHS - MODIFY THESE ACCORDING TO YOUR DATA STRUCTURE
# ===================================================================

# Base directories
DATA_DIR="raw_data"                    # Main data directory
GENOME_REF_DIR="${DATA_DIR}/genome_ref"  # Genome reference directory
FASTQ_DIR="${DATA_DIR}/fastq"      # Directory containing FASTQ files
BARCODE_DIR="${DATA_DIR}/barcodes" # Directory containing barcode files

# Reference Genome and Annotation Files
GENOME_FASTA="${GENOME_REF_DIR}/GCA_030064105.1_ASM3006410v1_genomic.fna"  # 
GFF3_FILE="${GENOME_REF_DIR}/genomic_annotation_R401.gff"                  # Path to annotation file

# Sequencing Read Files (Paired-End)
READ1="${FASTQ_DIR}/microSPLIT-600cells_S1_L001_R1_001.fastq.gz"  # Path to R1 file (contains barcodes)
READ2="${FASTQ_DIR}/microSPLIT-600cells_S1_L001_R2_001.fastq.gz"  # Path to R2 file (contains cDNA)

# Output Directory Configuration
GENOME_DIR="genome_index"           # Directory for STAR genome index
OUTPUT_DIR="starsolo_output"        # Directory for STARsolo output

# Barcode Files for Single-Cell Sequencing
# Specify barcode files for each round of barcoding
BARCODE_RD1="${BARCODE_DIR}/barcode_round1.txt"  # Barcodes for round 1
BARCODE_RD2="${BARCODE_DIR}/barcode_round2.txt"  # Barcodes for round 2
BARCODE_RD3="${BARCODE_DIR}/barcode_round3.txt"  # Barcodes for round 3
#BARCODE_R2R3="${BARCODE_DIR}/Barcode2-3.txt"    # Commented out alternative barcode file

# ===================================================================
# END OF PATH CONFIGURATION
# ===================================================================

# Create Necessary Output Directories
mkdir -p $GENOME_DIR $OUTPUT_DIR

# Input File Validation
# Check if all required input files exist before processing
echo "$(date) - Validating input files..."
for file in "$GENOME_FASTA" "$GFF3_FILE" "$READ1" "$READ2" "$BARCODE_RD1" "$BARCODE_RD2" "$BARCODE_RD3"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Missing input file: $file"
        exit 1
    fi
done

# Analysis Workflow
echo "$(date) - Starting STARsolo analysis" 

# 1. Convert GFF3 to GTF Format
# Use gffread to transform genome annotation from GFF3 to GTF
echo "$(date) - Converting GFF3 to GTF..."
GTF_FILE="${GFF3_FILE%.gff}.gtf"
gffread $GFF3_FILE -T -v -o $GTF_FILE
check_success "GFF3 to GTF conversion"

# 2. GTF File Correction
# Modify GTF to ensure compatibility with STAR
echo "$(date) - Correcting GTF file..."
GTF_FIXED="${GTF_FILE%.gtf}_fixed.gtf"
awk 'BEGIN{OFS=FS="\t"} {if ($3 == "CDS" || $3 == "transcript") $3 = "exon"; print}' $GTF_FILE > $GTF_FIXED
check_success "GTF file correction"

# 3. Genome Index Generation
# Create genome index for faster alignment
echo "$(date) - Generating genome index..."
STAR --runMode genomeGenerate \
    --runThreadN $THREADS \
    --genomeDir $GENOME_DIR \
    --genomeSAindexNbases 10 \
    --genomeFastaFiles $GENOME_FASTA \
    --sjdbGTFfile $GTF_FIXED
check_success "Genome index generation"

# 4. STARsolo Alignment and Quantification
# Perform single-cell RNA-seq alignment and gene counting with barcodes
# Cell barcode positions: 0_10_0_17 0_48_0_55 0_78_0_85 (3 rounds)
# UMI position in read: 0_0_0_9
echo "$(date) - Executing STARsolo alignment..."
STAR --genomeDir $GENOME_DIR \
    --runThreadN $THREADS \
    --readFilesIn $READ1 $READ2 \
    --readFilesCommand gunzip -c \
    --soloType CB_UMI_Complex \
    --soloCBposition 0_10_0_17 0_48_0_55 0_78_0_85 \
    --soloUMIposition 0_0_0_9 \
    --soloCBwhitelist $BARCODE_RD3 $BARCODE_RD2 $BARCODE_RD1 \
    --soloCBmatchWLtype 1MM \
    --soloUMIdedup 1MM_All \
    --soloFeatures Gene GeneFull \
    --soloMultiMappers Uniform \
    --outFilterScoreMinOverLread 0 \
    --outFilterMatchNminOverLread 0 \
    --outFilterMatchNmin 50 \
    --alignSJDBoverhangMin 1000 \
    --alignSJoverhangMin 1000 \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix "$OUTPUT_DIR/"
check_success "STARsolo alignment and quantification"

# 5. Output Validation
# Check if critical output files were generated
echo "$(date) - Verifying output files..."
if [ -f "$OUTPUT_DIR/Log.final.out" ] && [ -f "$OUTPUT_DIR/Solo.out/GeneFull/Raw/UniqueAndMult-Uniform.mtx" ]; then
    echo "Analysis completed successfully!"
    echo "Key files for downstream analysis:"
    echo "- $OUTPUT_DIR/Solo.out/GeneFull/barcodes.tsv"
    echo "- $OUTPUT_DIR/Solo.out/GeneFull/features.tsv"
    echo "- $OUTPUT_DIR/Solo.out/GeneFull/Raw/UniqueAndMult-Uniform.mtx"
    
    # Display alignment statistics
    echo "Alignment statistics:"
    grep "Uniquely mapped reads %" "$OUTPUT_DIR/Log.final.out"
    grep "Number of reads mapped to multiple loci" "$OUTPUT_DIR/Log.final.out"
    grep "Number of detected cells" "$OUTPUT_DIR/Solo.out/GeneFull/Summary.csv" 2>/dev/null || echo "Cell information not available"
else
    echo "ERROR: Some output files are missing!"
    exit 1
fi

# Create metadata file for traceability
echo "$(date) - Creating metadata file..."
{
    echo "# STARsolo analysis metadata"
    echo "Date: $(date)"
    echo "Genome: $GENOME_FASTA"
    echo "Annotation: $GFF3_FILE"
    echo "Read1: $READ1"
    echo "Read2: $READ2"
    echo "Command version: $(STAR --version)"
} > "$OUTPUT_DIR/analysis_metadata.txt"

echo "$(date) - Processing completed"

# ===================================================================   
# STARsolo Parameters Reference
# ===================================================================

# Barcode Positioning (--soloCBposition)
# 0_10_0_17 0_48_0_55 0_78_0_85: This indicates multiple barcode segments in specific read positions
# Each set of numbers represents start_length for different barcode rounds

# UMI Positioning (--soloUMIposition)
# 0_0_0_9: Defines the UMI location, here a 9-base sequence at the start of the read

# Whitelist and Matching (--soloCBwhitelist and --soloCBmatchWLtype)
# Uses predefined barcode lists from three rounds
# 1MM allows one mismatch for flexible barcode matching

# Deduplication (--soloUMIdedup)
# 1MM_All: Removes duplicates with one mismatch across all mapped locations
