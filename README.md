# GWAS Harmonization Pipeline

**Purpose:** Standardize heterogeneous GWAS summary statistics into a unified format with dual-build genomic coordinates.

---

## Overview

This **Nextflow DSL2 pipeline** harmonizes GWAS files from diverse sources (public consortia, biobanks) into a consistent column structure with both hg19 and hg38 coordinates. It handles:

- **Column name mapping** - Converts column names to standard format
- **Dual-build liftover** - Generates both original and lifted genome builds using MungeSumstats
- **Build merging** - Combines hg19/hg38 coordinates into single unified output
- **Quality filtering** - Removes sex chromosomes (X/Y/MT), non-biallelic variants, and duplicates
- **Allele harmonization** - Performs strand flipping and allele consistency checks
- **Comprehensive logging** - Records all transformations and filtering statistics per build

---

## Quick Start

### 1. Add Your GWAS File

```bash
# For non-biobank GWAS (for example)
mkdir -p NonBiobanks/raw_data/my_trait/
cp /path/to/original_gwas.txt.gz NonBiobanks/raw_data/my_trait/

# Add any metadata files to NonBiobanks/raw_data/my_trait/aux_files
```

### 2. Create JSON Configuration

```bash
# Copy template
cp /qc_scripts/json_configs/template.json /qc_scripts/json_configs/my_trait.json

# Edit with your GWAS column names
vi /qc_scripts/json_configs/my_trait.json
```

**Example JSON:**
```json
[
  {
    "trait": "trait_name",      # Name of the trait being analyzed
    "source": "NonBiobank",          # Does the GWAS come from a biobank or non-biobank source (value = NonBiobank or Biobank)
    "gwas_file": "/path/to/gwas/summary_statistics.txt",        # path to the GWAS summary statistics file
    "build": "37",        # genome build of the GWAS summary statistics (value = 37 or 38)
    "rsid_col": "SNP",      # rsID column (optional)
    "chr_col": "CHR",       # chromosome column
    "pos_col": "BP",    # position column
    "beta_col": "BETA",         # effect size column
    "se_col": "SE",     # standard error column
    "pval_col": "P",        # p-value column
    "A1": "A1",    # effect allele
    "A2": "A2",     # non-effect allele
    "eaf_col": "FRQ",     # effect allele frequency column
    "n_col": "N", # sample size column (optional)
    "info_col": "INFO"      # imputation quality column (optional)
  }
]
```

### 3. Run Harmonization

```bash
ml R
ml nextflow
cd qc_scripts/
nextflow run harmonize.nf \
  --json_file json_configs/my_trait.json \
  --output_dir /path/to/output \
  -c harmonize.config
```

**Output files:**
- `my_trait/my_trait_hg19.txt.gz` - hg19 coordinates
- `my_trait/my_trait_hg38.txt.gz` - hg38 coordinates
- `my_trait/my_trait_munged.txt.gz` - Merged dual-build file
- `my_trait/logs/original/` - Harmonization logs for original build
- `my_trait/logs/lifted/` - Harmonization logs for lifted build

### 4. Single-Build Mode (Optional)

```bash
# Skip build merging, output only hg19 and hg38 files separately
nextflow run harmonize.nf \
  --json_file json_configs/my_trait.json \
  --output_dir /path/to/output \
  --dual_build false \
  -c harmonize.config
```

---

## Directory Structure

```
GWASs/
├── NonBiobanks/                     
│   ├── raw_data/
│   │   ├── ad/                       # One directory per trait
│   │   │   ├── original.txt.gz       # Original GWAS file
│   │   │   
│   │   ├── bmi/
│   │   ├── cad/
│   │   └── ibd/
│   ├── qced_data/                    # Pipeline outputs
│   │   ├── ad/
│   │   │   ├── ad_hg19.txt.gz
│   │   │   ├── ad_hg38.txt.gz
│   │   │   ├── ad_munged.txt.gz      # Dual-build merged file
│   │   │   └── logs/
│   │   │       ├── original/         # MungeSumstats logs
│   │   │       └── lifted/           # MungeSumstats logs
│   │   └── ...
│   └── json_configs/                 # Configuration files
│       ├── ad.json
│       ├── bmi.json
│       └── ...
│
└── scripts/
    ├── harmonize.nf                  # Main Nextflow pipeline
    ├── harmonize.config              # Nextflow configuration
    ├── harmonize_gwas.R              # MungeSumstats wrapper
    ├── merge_builds.R                # Merge builds
    └── README.md                     # This file
```

---

## JSON Configuration Reference

### Required Fields

| Field | Description | Example |
|-------|-------------|---------|
| `trait` | Short trait identifier (used in output filenames) | `"ad"`, `"bmi"`, `"cad"` |
| `source` | Data source category | `"NonBiobank"`, `"Biobank_UKB", "Biobank_AoU"` |
| `gwas_file` | Path to raw GWAS file | `"../raw_data/ad/original.txt.gz"` |
| `build` | Original genome build | `"37"` or `"38"` |
| `chr_col` | Chromosome column
| `pos_col` | Base-pair position column 
| `A1` | Effect allele column
| `A2` | Non-effect allele column
| `beta_col` | Effect size column 
| `se_col` | Standard error column
| `pval_col` | P-value column

### Optional Fields

| Field | Description | Notes |
|-------|-------------|-------|
| `rsid_col` | rsID column | Will be included if present; variants identified by chr:pos:A1:A2 if missing |
| `eaf_col` | Effect allele frequency
| `n_col` | Sample size 
| `info_col` | Imputation quality score

> **Note:** Optional columns will be included in output if present in input file. By default, the munged sumstats file will include above columns, followed by any additional columns (appended to the end of the file).

---

## Standardized Output Format

### Merged file with both genome builds (`*_munged.txt.gz`)

| Column | Type | Description |
|--------|------|-------------|
| `SNP` | character | Variant identifier (rsID) |
| `CHR_hg19` | character | GRCh37 chromosome |
| `BP_hg19` | numeric | GRCh37 position |
| `CHR_hg38` | character | GRCh38 chromosome |
| `BP_hg38` | numeric | GRCh38 position |
| `A1` | character | Effect allele |
| `A2` | character | Non-effect allele |
| `BETA` | numeric | Effect size |
| `SE` | numeric | Standard error of BETA |
| `P` | numeric | P-value |
| `FRQ` | numeric | Effect allele frequency (optional) |
| `N` | numeric | Sample size (optional) |
| `INFO` | numeric | Imputation quality (optional) |

**File format:** Tab-delimited, gzip-compressed text (.txt.gz)

### Quality Filters Applied

1. Remove sex chromosomes (X, Y, MT)
2. Keep only biallelic SNPs
3. Remove duplicate variants (keeps first occurrence)
4. Perform allele flipping and strand alignment
5. Convert small/large/negative p-values to valid range

**Allele convention:** A1 = effect allele, A2 = non-effect allele  
*(Note: internally swapped during MungeSumstats processing to align with MungeSumstats requirements)*)*



---

## Pipeline Configuration

### Resource Requirements

Default settings in `harmonize.config`:

```nextflow
process {
    withName: HARMONIZE {
        cpus = 1
        memory = '8 GB'
        time = '2h'
    }
    
    withName: MERGE_BUILDS {
        cpus = 1
        memory = '4 GB'
        time = '30m'
    }
}
```

---

## FAQ

### Q: Can I skip the dual-build merging?

**A:** Yes, use `--dual_build false` to output only separate hg19/hg38 files.

### Q: What if my GWAS uses OR instead of BETA?

**A:** MungeSumstats automatically detects OR columns and converts to log(OR). Alternatively, pre-convert:

```bash
zcat original.txt.gz | \
    awk 'BEGIN {OFS="\t"} 
         NR==1 {print $0, "BETA"} 
         NR>1 {print $0, log($OR_col)}' | \
    gzip > with_beta.txt.gz
````

### Q: How are sex chromosomes handled?

**A:** Removed by default (rmv_chr = c("X", "Y", "MT")) for autosomal analyses. Modify harmonize_gwas.R to retain them.

### Q: Can I process build 38 GWASs?
**A:** Yes, set "build": "38" in JSON. Pipeline will lift to hg19 and merge both builds.



---

## Citation

If using this pipeline, please cite:

Murphy AE, et al. (2021). MungeSumstats: a Bioconductor package for the standardization and quality control of many GWAS summary statistics. *Bioinformatics*. doi: 10.1093/bioinformatics/btab665

---
