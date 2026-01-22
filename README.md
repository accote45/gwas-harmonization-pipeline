# GWAS Harmonization Pipeline

Nextflow pipeline for harmonizing GWAS summary statistics using MungeSumstats with dual-build support (hg19/hg38).

## Features

- Harmonizes GWAS summary statistics to standardized format
- Supports both build 37 (hg19) and build 38 (hg38)
- Automatic liftover between genome builds
- Merges dual builds into single output file
- Handles GWAS files with or without rsIDs
- Quality control and allele flipping

## Requirements

- Nextflow >= 21.04
- R >= 4.0
- R packages: MungeSumstats, jsonlite, data.table

## Quick Start

```bash
nextflow run harmonize.nf \
  --json_file json_configs/cad_nonbiobank.json \
  --output_dir /path/to/output