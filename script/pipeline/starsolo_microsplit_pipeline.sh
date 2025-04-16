#!/bin/bash
#SBATCH --job-name=STARsolo_Pipeline
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valentin.goupille@univ-rennes.fr
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=starsolo_pipeline_%j.out
#SBATCH --error=starsolo_pipeline_%j.err

# === STARsolo Pipeline Master Script ===
# This script runs the complete STARsolo analysis pipeline for microSPLIT scRNA-seq data

# Load configuration
echo "$(date) - Loading configuration..."

# Search for config.sh in multiple possible locations without absolute paths
CONFIG_FOUND=0
CONFIG_LOCATIONS=(
    "./config.sh"
    "script/pipeline/config.sh"
    "$(dirname "$0")/config.sh"
    "../config.sh"
    "../../config.sh"
)

for CONFIG_PATH in "${CONFIG_LOCATIONS[@]}"; do
    if [ -f "$CONFIG_PATH" ]; then
        echo "Found configuration file at: $CONFIG_PATH"
        source "$CONFIG_PATH"
        CONFIG_FOUND=1
        break
    fi
done

if [ $CONFIG_FOUND -eq 1 ]; then
    echo "Configuration loaded successfully."
else
    echo "ERROR: config.sh not found in any of the expected locations!"
    echo "Searched in:"
    for loc in "${CONFIG_LOCATIONS[@]}"; do
        echo "  - $loc"
    done
    echo "Please create the configuration file at one of these locations or create a symbolic link."
    exit 1
fi

# Verify BASE_PATH exists
if [ ! -d "$BASE_PATH" ]; then
    echo "ERROR: BASE_PATH directory does not exist: $BASE_PATH"
    echo "Please check the BASE_PATH variable in your config.sh file"
    exit 1
fi

# Create BASE_DIR early to avoid issues with conda environment
echo "$(date) - Creating base directory structure..."
mkdir -p "$BASE_DIR"
echo "Base directory created or already exists: $BASE_DIR"

# Note: SLURM parameters are now handled by the submit_pipeline.sh wrapper
# No need to update SLURM parameters here as they are passed via sbatch command
echo "Using SLURM parameters provided at submission time."

# Function to check if a script completed successfully
check_script() {
    if [ $? -ne 0 ]; then
        echo "ERROR: $1 failed! Pipeline stopped."
        exit 1
    else
        echo "$1 completed successfully."
    fi
}

# Step 0: Source Conda environment script
echo "$(date) - Setting up Conda environment..."
if [ -f "$CONDA_INIT_SCRIPT" ]; then
    echo "Sourcing Conda initialization script from: $CONDA_INIT_SCRIPT"
    . "$CONDA_INIT_SCRIPT"
else
    echo "ERROR: Conda initialization script not found at $CONDA_INIT_SCRIPT"
    echo "Please check the CONDA_INIT_SCRIPT variable in your config.sh file"
    exit 1
fi

# Step 1: Create STARsolo environment (if not already created)
echo "$(date) - Step 1: Creating STARsolo environment (if needed)..."

# First check if the parent directory of the conda environment exists
CONDA_PARENT_DIR=$(dirname "$CONDA_ENV_PATH")
if [ ! -d "$CONDA_PARENT_DIR" ]; then
    echo "Creating parent directory for conda environment: $CONDA_PARENT_DIR"
    mkdir -p "$CONDA_PARENT_DIR" || { echo "ERROR: Failed to create conda environment parent directory"; exit 1; }
fi

# Check if the environment already exists
if conda env list | grep -q "${CONDA_ENV_PATH##*/}"; then
    echo "STARsolo environment already exists. Trying to activate..."
    if ! conda activate "$CONDA_ENV_PATH"; then
        echo "ERROR: Cannot activate existing conda environment. It may be corrupted."
        echo "Removing and recreating the environment..."
        conda env remove -p "$CONDA_ENV_PATH" || echo "Warning: Failed to remove existing environment"
        CREATE_NEW_ENV=1
    else
        echo "Successfully activated existing environment."
        CREATE_NEW_ENV=0
    fi
else
    echo "STARsolo environment does not exist. Will create new environment."
    CREATE_NEW_ENV=1
fi

# Create a new environment if needed
if [ $CREATE_NEW_ENV -eq 1 ]; then
    echo "Creating STARsolo environment at: $CONDA_ENV_PATH"
    conda create -p "$CONDA_ENV_PATH" python=3.6 -y || { echo "ERROR: Failed to create conda environment"; exit 1; }
    
    echo "Installing required packages..."
    conda activate "$CONDA_ENV_PATH" || { echo "ERROR: Failed to activate the new conda environment"; exit 1; }
    conda install -c bioconda star=2.7.9 cufflinks=2.2.1 gffread=0.12.1 -y || { echo "ERROR: Failed to install required packages"; exit 1; }
    
    # Check installed tool versions
    echo "Checking installed tools versions:"
    STAR --version || { echo "ERROR: STAR not properly installed"; exit 1; }
    
    # Cufflinks doesn't support --version, just check if it's executable
    echo "Checking Cufflinks installation:"
    if ! command -v cufflinks &> /dev/null; then
        echo "ERROR: Cufflinks not found"
        exit 1
    else
        echo "Cufflinks is installed"
        # Display the first line of help output to confirm it works
        cufflinks 2>&1 | head -n 1
    fi
    
    # Check if gffread is available
    echo "Checking gffread installation:"
    if ! command -v gffread &> /dev/null; then
        echo "ERROR: gffread not found"
        exit 1
    else
        echo "gffread is installed"
    fi
fi

# Always activate the environment before continuing
echo "Activating conda environment..."
conda activate "$CONDA_ENV_PATH" || { echo "ERROR: Failed to activate conda environment at $CONDA_ENV_PATH"; exit 1; }

# Check that required commands are available
echo "Verifying required tools availability..."

# First check STAR and gffread which accept version options
for cmd in STAR gffread; do
    if ! command -v $cmd &> /dev/null; then
        echo "ERROR: $cmd command not found in the activated environment"
        echo "Please check your conda environment installation"
        exit 1
    fi
done

# Special check for cufflinks
if ! command -v cufflinks &> /dev/null; then
    echo "ERROR: cufflinks command not found in the activated environment"
    echo "Please check your conda environment installation"
    exit 1
else
    echo "Cufflinks is available"
fi

echo "Conda environment setup completed successfully."

# Step 2: Create directory structure and symbolic links
echo "$(date) - Step 2: Setting up directory structure and links..."

# Create directory structure
mkdir -p "$FASTQ_DIR" "$BARCODE_DIR" "$GENOME_REF_DIR" "$GENOME_ANNOTATION_DIR"

# Create symbolic links for FASTQ files
echo "Creating symbolic links for FASTQ files..."
ln -sf "$(realpath "$SOURCE_FASTQ/$SOURCE_FASTQ_R1")" "$FASTQ_DIR/$TARGET_FASTQ_R1"
ln -sf "$(realpath "$SOURCE_FASTQ/$SOURCE_FASTQ_R2")" "$FASTQ_DIR/$TARGET_FASTQ_R2"

# Create symbolic links for barcode files
echo "Creating symbolic links for barcode files..."
ln -sf "$(realpath "$SOURCE_BARCODES/barcode_round1.txt")" "$BARCODE_DIR/barcode_round1.txt"
ln -sf "$(realpath "$SOURCE_BARCODES/barcode_round2.txt")" "$BARCODE_DIR/barcode_round2.txt"
ln -sf "$(realpath "$SOURCE_BARCODES/barcode_round3.txt")" "$BARCODE_DIR/barcode_round3.txt"

# Create symbolic links for genome reference files
echo "Creating symbolic links for genome reference files..."
ln -sf "$(realpath "$SOURCE_GENOME_REF/$SOURCE_GENOME_FASTA")" "$GENOME_REF_DIR/$TARGET_GENOME_FASTA"

# Create symbolic links for genome annotation files
echo "Creating symbolic links for genome annotation files..."
ln -sf "$(realpath "$SOURCE_GENOME_ANNOTATION/$SOURCE_GENOME_GFF")" "$GENOME_ANNOTATION_DIR/$TARGET_GENOME_GFF"

check_script "Directory and link setup"

# Step 3: Run STARsolo analysis
echo "$(date) - Step 3: Running STARsolo analysis..."

# Create output directories
mkdir -p "$GENERAL_OUTPUT_DIR" "$GENOME_DIR" "$OUTPUT_STARSOLO_DIR"

# Activate STARsolo environment
conda activate "$CONDA_ENV_PATH"

# Function to check command success
check_success() {
    if [ $? -ne 0 ]; then
        echo "ERROR: Command failed - $1"
        exit 1
    fi
}

# Convert GFF3 to GTF Format
echo "$(date) - Converting GFF3 to GTF..."

# Ensure conda environment is still activated
if ! command -v gffread &> /dev/null; then
    echo "WARNING: gffread command not found. Trying to reactivate conda environment..."
    conda activate "$CONDA_ENV_PATH" || { 
        echo "ERROR: Failed to reactivate conda environment and gffread is not available."
        echo "Please ensure that gffread is installed in your conda environment."
        exit 1
    }
    
    # Check again after reactivation
    if ! command -v gffread &> /dev/null; then
        echo "ERROR: gffread command still not available after reactivating conda environment."
        echo "Please install gffread: conda install -c bioconda gffread"
        exit 1
    fi
fi

# Path for the output GTF file
GTF_FILE="${GENOME_ANNOTATION_DIR}/${TARGET_GENOME_GFF%.gff}.gtf"
echo "Converting GFF to GTF: ${GENOME_ANNOTATION_DIR}/${TARGET_GENOME_GFF} -> $GTF_FILE"

# Run conversion with -h option if -v doesn't work
gffread "${GENOME_ANNOTATION_DIR}/${TARGET_GENOME_GFF}" -T -o "$GTF_FILE"
if [ $? -ne 0 ]; then
    echo "ERROR: GFF3 to GTF conversion failed"
    exit 1
else
    echo "GFF3 to GTF conversion completed successfully"
fi

# Check that the GTF file was created
if [ ! -f "$GTF_FILE" ]; then
    echo "ERROR: GTF file was not created: $GTF_FILE"
    exit 1
fi

# GTF File Correction
echo "$(date) - Correcting GTF file..."
GTF_FIXED="${GTF_FILE%.gtf}_fixed.gtf"
awk 'BEGIN{OFS=FS="\t"} {if ($3 == "CDS" || $3 == "transcript") $3 = "exon"; print}' "$GTF_FILE" > "$GTF_FIXED"
check_success "GTF file correction"

# Genome Index Generation
echo "$(date) - Generating genome index..."
STAR --runMode genomeGenerate \
    --runThreadN "$THREADS" \
    --genomeDir "$GENOME_DIR" \
    --genomeSAindexNbases 10 \
    --genomeFastaFiles "${GENOME_REF_DIR}/${TARGET_GENOME_FASTA}" \
    --sjdbGTFfile "$GTF_FIXED"
check_success "Genome index generation"

# STARsolo Alignment and Quantification
echo "$(date) - Executing STARsolo alignment..."
STAR --genomeDir "$GENOME_DIR" \
    --runThreadN "$THREADS" \
    --readFilesIn "${FASTQ_DIR}/${TARGET_FASTQ_R1}" "${FASTQ_DIR}/${TARGET_FASTQ_R2}" \
    --readFilesCommand gunzip -c \
    --soloType CB_UMI_Complex \
    --soloCBposition 0_10_0_17 0_48_0_55 0_78_0_85 \
    --soloUMIposition 0_0_0_9 \
    --soloCBwhitelist "${BARCODE_DIR}/barcode_round3.txt" "${BARCODE_DIR}/barcode_round2.txt" "${BARCODE_DIR}/barcode_round1.txt" \
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
    --outFileNamePrefix "${OUTPUT_STARSOLO_DIR}/"
check_success "STARsolo alignment and quantification"

# Output Validation
echo "$(date) - Verifying output files..."

# Check for log file
if [ ! -f "${OUTPUT_STARSOLO_DIR}/Log.final.out" ]; then
    echo "ERROR: Log file not found!"
    exit 1
fi

# Check for matrix files in both Gene and GeneFull features
GENE_MTX="${OUTPUT_STARSOLO_DIR}/Solo.out/Gene/raw/UniqueAndMult-Uniform.mtx"
GENEFULL_MTX="${OUTPUT_STARSOLO_DIR}/Solo.out/GeneFull/raw/UniqueAndMult-Uniform.mtx"
BARCODE_FILE="${OUTPUT_STARSOLO_DIR}/Solo.out/GeneFull/raw/barcodes.tsv"
FEATURE_FILE="${OUTPUT_STARSOLO_DIR}/Solo.out/GeneFull/raw/features.tsv"

MISSING_FILES=0
for file in "$GENE_MTX" "$GENEFULL_MTX" "$BARCODE_FILE" "$FEATURE_FILE"; do
    if [ ! -f "$file" ]; then
        echo "WARNING: Missing output file: $file"
        MISSING_FILES=$((MISSING_FILES+1))
    fi
done

if [ $MISSING_FILES -gt 0 ]; then
    echo "WARNING: $MISSING_FILES output files are missing!"
else
    echo "Analysis completed successfully!"
    echo "Key files for downstream analysis:"
    echo "- $BARCODE_FILE"
    echo "- $FEATURE_FILE"
    echo "- $GENEFULL_MTX"
    echo "- $GENE_MTX"
    
    # Display detailed alignment statistics
    echo -e "\nAlignment statistics:"
    echo "-----------------------"
    echo "$(grep "Number of input reads" "${OUTPUT_STARSOLO_DIR}/Log.final.out")"
    echo "$(grep "Average input read length" "${OUTPUT_STARSOLO_DIR}/Log.final.out")"
    echo "$(grep "Uniquely mapped reads number" "${OUTPUT_STARSOLO_DIR}/Log.final.out")"
    echo "$(grep "Uniquely mapped reads %" "${OUTPUT_STARSOLO_DIR}/Log.final.out")"
    echo "$(grep "Number of reads mapped to multiple loci" "${OUTPUT_STARSOLO_DIR}/Log.final.out")"
    echo "$(grep "% of reads mapped to multiple loci" "${OUTPUT_STARSOLO_DIR}/Log.final.out")"
    echo "$(grep "Number of reads unmapped: other" "${OUTPUT_STARSOLO_DIR}/Log.final.out")"
    echo "$(grep "% of reads unmapped: other" "${OUTPUT_STARSOLO_DIR}/Log.final.out")"
    echo "$(grep "Mismatch rate per base, %" "${OUTPUT_STARSOLO_DIR}/Log.final.out")"
fi

echo "$(date) - Pipeline completed!" 