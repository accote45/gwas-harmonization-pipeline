#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript merge_builds.R <trait> <hg19_file> <hg38_file> <output_file>")
}

trait <- args[1]
hg19_file <- args[2]
hg38_file <- args[3]
output_file <- args[4]

cat("Merging builds for", trait, "...\n")

hg19 <- fread(hg19_file)
hg38 <- fread(hg38_file)

# Swap A1 and A2 to convert from MungeSumstats convention (A2=effect) 
# to desired convention (A1=effect)
if ("A1" %in% names(hg19) && "A2" %in% names(hg19)) {
  setnames(hg19, c("A1", "A2"), c("A2_temp", "A1_temp"))
  setnames(hg19, c("A1_temp", "A2_temp"), c("A1", "A2"))
}

if ("A1" %in% names(hg38) && "A2" %in% names(hg38)) {
  setnames(hg38, c("A1", "A2"), c("A2_temp", "A1_temp"))
  setnames(hg38, c("A1_temp", "A2_temp"), c("A1", "A2"))
}

setnames(hg19, c("CHR", "BP"), c("CHR_hg19", "BP_hg19"))
setnames(hg38, c("CHR", "BP"), c("CHR_hg38", "BP_hg38"))

if ("SNP" %in% names(hg19)) {
  merge_key <- "SNP"
} else {
  hg19[, merge_key := paste(A1, A2, CHR_hg19, BP_hg19, sep = "_")]
  hg38[, merge_key := paste(A1, A2, CHR_hg38, BP_hg38, sep = "_")]
  merge_key <- "merge_key"
}

hg38_coords <- hg38[, c(merge_key, "CHR_hg38", "BP_hg38"), with = FALSE]
merged <- merge(hg19, hg38_coords, by = merge_key, all.x = TRUE)

if (merge_key == "merge_key") {
  merged[, merge_key := NULL]
}

# Define proper column ordering with A1 before A2
coord_cols <- c("SNP", "CHR_hg19", "BP_hg19", "CHR_hg38", "BP_hg38", "A1", "A2")
coord_cols <- coord_cols[coord_cols %in% names(merged)]
other_cols <- setdiff(names(merged), coord_cols)
setcolorder(merged, c(coord_cols, other_cols))

fwrite(merged, output_file, sep = "\t")

cat("✓ Merged", nrow(merged), "variants\n")
cat("  Output:", output_file, "\n")