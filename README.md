# STARsolo Pipeline pour l'Analyse de Données microSPLIT scRNA-seq

Pipeline pour l'analyse des données de séquençage de cellules uniques (scRNA-seq) avec la technologie microSPLIT utilisant STARsolo.

## Description

Ce projet contient une série de scripts pour analyser des données de RNA-seq unicellulaire générées avec la technologie microSPLIT. Il utilise STARsolo, une implémentation de la pipeline STAR dédiée aux données unicellulaires, pour aligner les reads et quantifier l'expression génique.

## Structure du Projet

```
.
├── create_env_starsolo.sh    # Script pour créer l'environnement conda avec STARsolo
├── create_links.sh           # Script pour préparer la structure des répertoires et créer des liens symboliques
├── run_starsolo_DOL_microsplit.sh  # Script principal d'analyse STARsolo
└── raw_data/                 # Répertoire contenant les données brutes
    ├── fastq/                # Fichiers FASTQ (reads)
    ├── barcodes/             # Fichiers de barcodes
    ├── genome_ref/           # Fichiers de référence du génome
    └── genome_annotation/    # Fichiers d'annotation du génome
```

## Prérequis

- Conda/Miniconda
- Accès à un cluster de calcul (le script est configuré pour SLURM)
- Données brutes organisées comme décrit dans la section suivante

## Installation

1. Créez l'environnement conda avec STARsolo:

```bash
sbatch create_env_starsolo.sh
```

2. Organisez vos données brutes dans la structure appropriée:

```
raw_data/
├── fastq/
│   ├── microSPLIT-600cells_S1_L001_R1_001.fastq.gz
│   └── microSPLIT-600cells_S1_L001_R2_001.fastq.gz
├── barcodes/
│   ├── barcode_round1.txt
│   ├── barcode_round2.txt
│   └── barcode_round3.txt
├── genome_ref/
│   └── GCA_030064105.1_ASM3006410v1_genomic.fna
└── genome_annotation/
    └── GCA_030064105.1_ASM3006410v1_genomic.gff
```

3. Créez les liens symboliques nécessaires:

```bash
sbatch create_links.sh
```

## Utilisation

Lancer l'analyse STARsolo:

```bash
sbatch run_starsolo_DOL_microsplit.sh
```

Ce script effectue les étapes suivantes:
1. Conversion du fichier GFF3 en GTF
2. Correction du fichier GTF pour compatibilité avec STAR
3. Génération de l'index du génome
4. Alignement et quantification avec STARsolo
5. Vérification et rapport des résultats

## Paramètres Personnalisables

Les paramètres principaux peuvent être modifiés dans le script `run_starsolo_DOL_microsplit.sh`:

- `THREADS`: Nombre de threads à utiliser
- `BASE_DIR`: Répertoire de base pour l'analyse
- Chemins vers les fichiers d'entrée et de sortie
- Positions des barcodes cellulaires et UMI dans les reads

## Résultats

Les résultats seront générés dans le répertoire `starsolo_script_DOL_microsplit/Output_DOL_microsplit_starsolo/` avec la structure suivante:

```
Output_DOL_microsplit_starsolo/
├── genome_index/          # Index du génome généré par STAR
└── starsolo_output/       # Résultats de l'analyse STARsolo
    ├── Log.final.out      # Résumé des statistiques d'alignement
    └── Solo.out/          # Résultats spécifiques à l'analyse unicellulaire
        ├── Gene/          # Matrices d'expression au niveau des gènes
        └── GeneFull/      # Matrices d'expression avec tous les reads
```

## Remarques

- L'analyse est configurée spécifiquement pour la technologie microSPLIT avec 3 rounds de barcoding cellulaire.
- Les paramètres d'alignement STAR sont optimisés pour des données microSPLIT. 