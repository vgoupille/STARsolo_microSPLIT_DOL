# STARsolo Pipeline for microSPLIT scRNA-seq Data Analysis

Pipeline for analyzing single-cell RNA sequencing (scRNA-seq) data generated with microSPLIT technology using STARsolo.

## Getting Started

### Access Cluster and Clone Repository

1. Connect to your computational cluster:
```bash
ssh your_username@your_cluster_address
```

2. Navigate to your preferred working directory:
```bash
cd /path/to/your/working/directory
# For example:
# cd /home/user/projects/scRNAseq_analysis
```

3. Clone the repository and set up the environment:
```bash
# Clone the repository
git clone https://github.com/vgoupille/STARsolo_microSPLIT_DOL.git
cd STARsolo_microSPLIT_DOL

# Make the submission script executable
chmod +x submit_pipeline.sh
```

## Description

This project contains a series of scripts to analyze single-cell RNA-seq data generated with microSPLIT technology. It uses STARsolo, an implementation of the STAR pipeline dedicated to single-cell data, to align reads and quantify gene expression. The pipeline generates count matrices that are essential for downstream analysis using tools like Seurat (R) or Scanpy (Python).

### About microSPLIT Technology

This pipeline is specifically designed for data generated using the microSPLIT protocol, as described in [Kuchina et al. (2021)](https://www.science.org/doi/10.1126/science.aba5257),
[Gaiser et al. (2024)](https://www.nature.com/articles/s41596-024-01007-w) which uses a split-pool barcoding approach for single-cell sequencing for bacteria.

### Generated Output Matrices

STARsolo generates several types of count matrices depending on different mapping and counting strategies:

- **Mapping strategies**: 
  - `Unique`: Only uniquely mapped reads
  - `UniqueAndMult-Uniform`: Both uniquely mapped reads and multimapped reads (with uniform weight distribution)

- **Feature types**:
  - `Gene`: Counts based on reads mapped to exons only
  - `GeneFull`: Counts including reads mapped to both exons and introns

- **Filtering options**:
  - `raw`: All barcodes, including potential background
  - `filtered`: Only cells passing quality thresholds

The current implementation uses `UniqueAndMult-Uniform` mapping with both `Gene` and `GeneFull` feature counting to maximize data recovery, especially for species with incomplete annotations.

### Key Parameters

The pipeline uses optimized parameters for microSPLIT data:
- Complex barcode setup with 3 rounds of barcoding
- CB (Cell Barcode) positions: 0_10_0_17 0_48_0_55 0_78_0_85 (corresponding to the 3 rounds)
- UMI position: 0_0_0_9
- 1 mismatch allowed in barcode matching
- Optimized alignment parameters for increased sensitivity

For more details on STARsolo parameters and options, see the [official documentation](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md).

## Project Structure

```
.
├── config.sh                      # Configuration file with all parameters
├── submit_pipeline.sh             # Wrapper script to submit the pipeline with custom parameters
├── starsolo_microsplit_pipeline.sh # Main script implementing the complete analysis pipeline
├── create_env_starsolo.sh         # Script to create conda environment with STARsolo (legacy)
├── create_links.sh                # Script to prepare directory structure and links (legacy)
├── run_starsolo_DOL_microsplit.sh # STARsolo analysis script (legacy)
└── raw_data/                      # Directory containing raw data
    ├── fastq/                     # FASTQ files (reads)
    ├── barcodes/                  # Barcode files
    ├── genome_ref/                # Genome reference files
    └── genome_annotation/         # Genome annotation files
```

## Prerequisites

- Conda/Miniconda
- Access to a computing cluster (script is configured for SLURM)
- Raw data organized as described in the configuration file

## Installation and Usage

### New Wrapper Approach (Recommended)

1. Edit the `config.sh` file to match your data paths and resource requirements:

```bash
nano config.sh
```

2. Run the pipeline with default parameters from config.sh:

```bash
./submit_pipeline.sh
```

3. Or customize SLURM parameters at submission time:

```bash
./submit_pipeline.sh --threads 32 --memory 64G --runtime 12:00:00
```

4. For all available options:

```bash
./submit_pipeline.sh --help
```

### Direct Submission (Alternative)

You can still submit the pipeline script directly with SLURM, but you won't be able to override parameters as easily:

```bash
sbatch starsolo_microsplit_pipeline.sh
```

### Legacy Approach (Individual Scripts)

If you prefer, you can still run each step individually:

1. Create the conda environment with STARsolo:

```bash
sbatch create_env_starsolo.sh
```

2. Create the necessary symbolic links:

```bash
sbatch create_links.sh
```

3. Launch the STARsolo analysis:

```bash
sbatch run_starsolo_DOL_microsplit.sh
```

## Customizable Parameters

### Configuration File (config.sh)

The configuration file (`config.sh`) contains all parameters needed to run the pipeline. **You must review and modify this file before running the pipeline** to ensure it matches your data and computing environment.

#### Essential Parameters (Must Configure)

- **Data Paths**:
  - `SOURCE_FASTQ`: Absolute path to the directory containing your FASTQ files
  - `SOURCE_BARCODES`: Absolute path to the directory containing barcode files
  - `SOURCE_GENOME_REF`: Absolute path to the directory containing genome reference FASTA
  - `SOURCE_GENOME_ANNOTATION`: Absolute path to the directory containing genome annotation (GFF)

- **Input File Names**:
  - `SOURCE_FASTQ_R1`: Filename of the R1 FASTQ file (reads)
  - `SOURCE_FASTQ_R2`: Filename of the R2 FASTQ file (reads)
  - `SOURCE_GENOME_FASTA`: Filename of the genome reference FASTA
  - `SOURCE_GENOME_GFF`: Filename of the genome annotation file (GFF format)

#### Recommended to Configure

- **Resources** (can also be overridden through submit_pipeline.sh):
  - `THREADS`: Number of CPU threads to use (default: 16, recommended: 32 or 64 for large datasets)
  - `MEMORY`: Memory allocation (default: 16G, recommended: 32G or 64G for large genomes)
  - `MAX_RUNTIME`: Maximum job runtime in format HH:MM:SS (default: 10:00:00)

- **Analysis Identification**:
  - `EMAIL`: Email address for job notifications
  - `ANALYSIS_NAME`: Name of the analysis (will be used for job name)
  - `BASE_DIR`: Directory where all output will be stored

#### Optional Parameters

- **Environment**:
  - `CONDA_ENV_PATH`: Path where the conda environment will be created

- **Output File Names**:
  - `TARGET_FASTQ_R1`: Internal symbolic link name for R1 FASTQ
  - `TARGET_FASTQ_R2`: Internal symbolic link name for R2 FASTQ
  - `TARGET_GENOME_FASTA`: Internal symbolic link name for genome FASTA
  - `TARGET_GENOME_GFF`: Internal symbolic link name for genome GFF

- **Internal Directory Structure**:
  - Several variables define the internal directory structure (typically don't need modification)

### Command-line Options (submit_pipeline.sh)

The following parameters can be overridden at submission time:

- `--threads N`: Number of CPU threads to use
- `--memory SIZE`: Memory allocation (e.g., 32G, 64G)
- `--runtime HH:MM:SS`: Maximum runtime
- `--email ADDRESS`: Notification email
- `--name JOBNAME`: SLURM job name

## Pipeline Steps

The pipeline performs the following steps:
1. Setting up the STARsolo environment
2. Creating the directory structure and symbolic links
3. Converting GFF3 file to GTF
4. Correcting GTF file for compatibility with STAR
5. Generating genome index
6. Alignment and quantification with STARsolo
7. Verification and reporting of results

## Results

The results will be generated in the configured output directory with the following structure:

```
Analysis_STARsolo_microsplit/
├── genome_index/          # Genome index generated by STAR
└── starsolo_output/       # STARsolo analysis results
    ├── Log.final.out      # Summary of alignment statistics
    └── Solo.out/          # Results specific to single-cell analysis
        ├── Gene/          # Expression matrices at gene level
        └── GeneFull/      # Expression matrices with all reads
```

## Notes

- The analysis is specifically configured for microSPLIT technology with 3 rounds of cell barcoding.
- STAR alignment parameters are optimized for microSPLIT data.
- The wrapper approach (submit_pipeline.sh) provides flexibility to customize resources without modifying the pipeline script.

## Author

**Valentin Goupille** - Université de Rennes

For questions or issues related to this pipeline, please open an issue on GitHub or contact the author directly. 