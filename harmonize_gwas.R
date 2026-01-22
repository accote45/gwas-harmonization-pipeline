#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(MungeSumstats)
  library(jsonlite)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript harmonize_gwas.R <json_file> <output_file> <build_type>")
}

json_file <- args[1]
output_file <- args[2]
build_type <- args[3]  # 'original' or 'lifted'

# Load config
cfg <- fromJSON(json_file)
if (is.list(cfg) && !"trait" %in% names(cfg)) cfg <- cfg[[1]]

original_build <- if (!is.null(cfg$build) && cfg$build == "38") "38" else "37"

# Determine if rsIDs are present
has_rsid <- !is.null(cfg$rsid_col) && nzchar(cfg$rsid_col)

# Column mapping
  # A1 and A2 are swapped to match MungeSumstats expectations
  # (MungeSumstats expects A1 = non-effect allele, A2 = effect allele)
col_map <- list(
  SNP = cfg$rsid_col, CHR = cfg$chr_col, BP = cfg$pos_col,
  A1 = cfg$A2, A2 = cfg$A1, BETA = cfg$beta_col,
  SE = cfg$se_col, P = cfg$pval_col, FRQ = cfg$eaf_col,
  N = cfg$n_col, INFO = cfg$info_col
)
col_map <- col_map[!sapply(col_map, is.null)]

column_mapping <- data.frame(
  Uncorrected = toupper(unname(unlist(col_map))),
  Corrected = names(col_map),
  stringsAsFactors = FALSE
)

# Determine genome builds
ref_genome <- if (original_build == "37") "GRCh37" else "GRCh38"
convert_genome <- if (build_type == "lifted") {
  if (original_build == "37") "GRCh38" else "GRCh37"
} else NULL

# Create log directory - simpler structure
log_dir <- "logs"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# Harmonize
format_sumstats(
  path = cfg$gwas_file,
  ref_genome = ref_genome,
  convert_ref_genome = convert_genome,
  
  INFO_filter = 0,
  FRQ_filter = 0,
  rmv_chr = c("X", "Y", "MT"),
  bi_allelic_filter = TRUE,
  
  allele_flip_check = TRUE,
  allele_flip_drop = TRUE,
  allele_flip_z = TRUE,
  allele_flip_frq = TRUE,
  
  check_dups = TRUE,
  sort_coordinates = TRUE,
  N_dropNA = FALSE,
  
  save_path = output_file,
  write_vcf = FALSE,
  
  log_folder_ind = FALSE,
  log_folder = log_dir,
  log_mungesumstats_msgs = TRUE,
  
  convert_small_p = TRUE,
  convert_large_p = TRUE,
  convert_neg_p = TRUE,
  
  snp_ids_are_rs_ids = has_rsid,
  
  mapping_file = column_mapping,
  force_new = TRUE
)

cat("✓ Completed", build_type, "harmonization\n")