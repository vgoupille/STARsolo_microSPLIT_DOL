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
if [ -f "config.sh" ]; then
    source config.sh
    echo "Configuration loaded successfully."
else
    echo "ERROR: config.sh not found! Please create the configuration file first."
    exit 1
fi

# Create BASE_DIR early to avoid issues with conda environment
echo "$(date) - Creating base directory structure..."
mkdir -p "$BASE_DIR"
echo "Base directory created or already exists: $BASE_DIR"

# Update SLURM parameters from config
# Note: These won't take effect for the current job, but are updated for consistency
sed -i "s/--job-name=.*/--job-name=${ANALYSIS_NAME}/g" "$0"
sed -i "s/--mail-user=.*/--mail-user=${EMAIL}/g" "$0"
sed -i "s/--cpus-per-task=.*/--cpus-per-task=${THREADS}/g" "$0"
sed -i "s/--mem=.*/--mem=${MEMORY}/g" "$0"
sed -i "s/--time=.*/--time=${MAX_RUNTIME}/g" "$0"

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
. /local/env/envconda.sh

# Step 1: Create STARsolo environment (if not already created)
echo "$(date) - Step 1: Creating STARsolo environment (if needed)..."
if conda env list | grep -q "${CONDA_ENV_PATH##*/}"; then
    echo "STARsolo environment already exists. Skipping creation."
else
    echo "Creating STARsolo environment..."
    # Extract commands from create_env_starsolo.sh and run them
    conda create -p "$CONDA_ENV_PATH" python=3.6 -y
    check_script "Environment creation"
    
    conda activate "$CONDA_ENV_PATH"
    conda install -c bioconda star=2.7.9 cufflinks=2.2.1 -y
    check_script "Package installation"
    
    # Print versions
    echo "STAR version:"
    STAR --version
    echo "Cufflinks version:"
    cufflinks --version
fi

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
GTF_FILE="${GENOME_ANNOTATION_DIR}/${TARGET_GENOME_GFF%.gff}.gtf"
gffread "${GENOME_ANNOTATION_DIR}/${TARGET_GENOME_GFF}" -T -v -o "$GTF_FILE"
check_success "GFF3 to GTF conversion"

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
BARCODE_FILE="${OUTPUT_STARSOLO_DIR}/Solo.out/GeneFull/barcodes.tsv"
FEATURE_FILE="${OUTPUT_STARSOLO_DIR}/Solo.out/GeneFull/features.tsv"

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