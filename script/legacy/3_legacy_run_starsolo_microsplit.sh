#!/bin/bash

# SLURM Job Configuration
#SBATCH --job-name=STARsolo_DOL_microSPLIT  # Job name
#SBATCH --mail-user=valentin.goupille@univ-rennes.fr
#SBATCH --mail-type=BEGIN,END,FAIL          # Email notifications
#SBATCH --ntasks=1                   # Number of tasks
#SBATCH --cpus-per-task=16  #16           # CPUs per task
#SBATCH --mem=16G  #64G                    # Memory allocation
#SBATCH --time=4:00:00              # Maximum job runtime
#SBATCH --output=starsolo_%j.out     # Standard output log file
#SBATCH --error=starsolo_%j.err      # Standard error log file

# Conda Environment Setup
# Source the Conda initialization script
. /local/env/envconda.sh

# Activate the specific Conda environment for STARsolo
conda activate 5_environnements/env_STARsolo  # Path to STARsolo environment

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
#choose the location of the folder for the Analysis : You can change it
BASE_DIR="Analysis_STARsolo_microsplit"
DATA_DIR="${BASE_DIR}/raw_data"                    # Main data directory
GENOME_REF_DIR="${DATA_DIR}/genome_ref"  # Genome reference directory
GENOME_ANNOTATION_DIR="${DATA_DIR}/genome_annotation" # Directory containing genome annotation
FASTQ_DIR="${DATA_DIR}/fastq"      # Directory containing FASTQ files
BARCODE_DIR="${DATA_DIR}/barcodes" # Directory containing barcode files


# Reference Genome and Annotation Files
GENOME_FASTA="${GENOME_REF_DIR}/genome_ref_PsR401.fna"  #path to genome fasta file (.fna)
GFF3_FILE="${GENOME_ANNOTATION_DIR}/genome_annotation_PsR401.gff"  #path to genome annotation file (.gff)

# Sequencing Read Files (Paired-End) (microsplit)
READ1="${FASTQ_DIR}/microSPLIT-R1.fastq.gz"  # Path to R1 file (contains barcodes)
READ2="${FASTQ_DIR}/microSPLIT-R2.fastq.gz"  # Path to R2 file (contains cDNA)

# Output Directory Configuration
GENERAL_OUTPUT_DIR="${BASE_DIR}/Output_DOL_microsplit_starsolo"
GENOME_DIR="${GENERAL_OUTPUT_DIR}/genome_index"           # Directory for STAR genome index
OUTPUT_STARSOLO_DIR="${GENERAL_OUTPUT_DIR}/starsolo_output"        # Directory for STARsolo output

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
mkdir -p $GENERAL_OUTPUT_DIR
mkdir -p $GENOME_DIR $OUTPUT_STARSOLO_DIR

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
    --outFileNamePrefix "$OUTPUT_STARSOLO_DIR/"
check_success "STARsolo alignment and quantification"

# 5. Output Validation
# Check if critical output files were generated and provide detailed metrics
echo "$(date) - Verifying output files..."

# Check for log file
if [ ! -f "$OUTPUT_STARSOLO_DIR/Log.final.out" ]; then
    echo "ERROR: Log file not found!"
    exit 1
fi

# Check for matrix files in both Gene and GeneFull features
GENE_MTX="$OUTPUT_STARSOLO_DIR/Solo.out/Gene/raw/UniqueAndMult-Uniform.mtx"
GENEFULL_MTX="$OUTPUT_STARSOLO_DIR/Solo.out/GeneFull/raw/UniqueAndMult-Uniform.mtx"
BARCODE_FILE="$OUTPUT_STARSOLO_DIR/Solo.out/GeneFull/barcodes.tsv"
FEATURE_FILE="$OUTPUT_STARSOLO_DIR/Solo.out/GeneFull/features.tsv"

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
    echo "$(grep "Number of input reads" "$OUTPUT_STARSOLO_DIR/Log.final.out")"
    echo "$(grep "Average input read length" "$OUTPUT_STARSOLO_DIR/Log.final.out")"
    echo "$(grep "Uniquely mapped reads number" "$OUTPUT_STARSOLO_DIR/Log.final.out")"
    echo "$(grep "Uniquely mapped reads %" "$OUTPUT_STARSOLO_DIR/Log.final.out")"
    echo "$(grep "Number of reads mapped to multiple loci" "$OUTPUT_STARSOLO_DIR/Log.final.out")"
    echo "$(grep "% of reads mapped to multiple loci" "$OUTPUT_STARSOLO_DIR/Log.final.out")"
    echo "$(grep "Number of reads unmapped: other" "$OUTPUT_STARSOLO_DIR/Log.final.out")"
    echo "$(grep "% of reads unmapped: other" "$OUTPUT_STARSOLO_DIR/Log.final.out")"
    echo "$(grep "Mismatch rate per base, %" "$OUTPUT_STARSOLO_DIR/Log.final.out")"
    
    # Display mapping speed
    echo -e "\nProcessing speed:"
    echo "--------------------"
    echo "$(grep "Mapping speed, Million of reads per hour" "$OUTPUT_STARSOLO_DIR/Log.final.out")"
    echo "$(grep "Started job on" "$OUTPUT_STARSOLO_DIR/Log.final.out")"
    echo "$(grep "Finished on" "$OUTPUT_STARSOLO_DIR/Log.final.out")"
    
    # Display cell statistics if available
    echo -e "\nCell statistics:"
    echo "----------------"
    if [ -f "$OUTPUT_STARSOLO_DIR/Solo.out/GeneFull/Summary.csv" ]; then
        SUMMARY_FILE="$OUTPUT_STARSOLO_DIR/Solo.out/GeneFull/Summary.csv"
        
        # Extract key metrics from Summary.csv
        echo "# Cell metrics from Summary.csv:"
        echo "--------------------------------"
        
        # Function to extract value from Summary.csv
        extract_metric() {
            grep "^$1," "$SUMMARY_FILE" | cut -d',' -f2
        }
        
        # Display important cell metrics
        echo "Estimated Number of Cells: $(extract_metric "Estimated Number of Cells")"
        echo "Reads With Valid Barcodes: $(extract_metric "Reads With Valid Barcodes") (fraction)"
        echo "Sequencing Saturation: $(extract_metric "Sequencing Saturation") (fraction)"
        echo "Reads Mapped to GeneFull: $(extract_metric "Reads Mapped to GeneFull: Unique+Multipe GeneFull") (fraction)"
        echo "Mean Reads per Cell: $(extract_metric "Mean Reads per Cell")"
        echo "Median Reads per Cell: $(extract_metric "Median Reads per Cell")"
        echo "UMIs in Cells: $(extract_metric "UMIs in Cells")"
        echo "Mean UMI per Cell: $(extract_metric "Mean UMI per Cell")"
        echo "Median UMI per Cell: $(extract_metric "Median UMI per Cell")"
        echo "Mean GeneFull per Cell: $(extract_metric "Mean GeneFull per Cell")"
        echo "Median GeneFull per Cell: $(extract_metric "Median GeneFull per Cell")"
        echo "Total GeneFull Detected: $(extract_metric "Total GeneFull Detected")"
        
        # Also show simple counts from files
        echo -e "\n# File-based counts:"
        echo "Barcodes count: $(wc -l < "$BARCODE_FILE")"
        echo "Features count: $(wc -l < "$FEATURE_FILE")"
    else
        echo "WARNING: Summary.csv file not found at $OUTPUT_STARSOLO_DIR/Solo.out/GeneFull/Summary.csv"
        echo "Estimated number of cells: $(wc -l < "$BARCODE_FILE")"
        echo "Number of genes/features: $(wc -l < "$FEATURE_FILE")"
    fi
fi

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
