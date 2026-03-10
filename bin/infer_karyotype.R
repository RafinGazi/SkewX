#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
vcf_file <- args[1]
individual <- args[2]

library(tidyverse)

# Read VCF directly - skip comment lines starting with ##
vcf_lines <- readLines(pipe(paste("zcat", vcf_file)))
header_line <- which(startsWith(vcf_lines, "#CHROM"))
vcf <- read_tsv(
    I(vcf_lines[header_line:length(vcf_lines)]),
    show_col_types=FALSE
)

# Filter chrX only
colnames(vcf)[1] <- "CHROM"
vcf <- vcf %>% filter(CHROM == "chrX")

# Extract genotype field (last column, before first colon)
gt_col <- vcf[[ncol(vcf)]]
gt <- sapply(strsplit(gt_col, ":"), `[`, 1)

# Count het and total calls
total <- sum(gt %in% c("0|0","0|1","1|0","1|1","0/0","0/1","1/0","1/1"))
het   <- sum(gt %in% c("0|1","1|0","0/1","1/0"))
het_ratio <- ifelse(total > 0, het / total, 0)

# Infer karyotype
karyotype <- case_when(
    het_ratio > 0.25 ~ "XX",
    het_ratio < 0.10 ~ "XY",
    TRUE             ~ "uncertain"
)
confidence <- case_when(
    het_ratio > 0.35 | het_ratio < 0.05 ~ "HIGH",
    het_ratio > 0.25 | het_ratio < 0.10 ~ "MEDIUM",
    TRUE                                 ~ "LOW"
)

# Write output
result <- tibble(
    individual  = individual,
    karyotype   = karyotype,
    het_ratio   = round(het_ratio, 4),
    het_count   = het,
    total_count = total,
    confidence  = confidence
)

write_tsv(result, paste0(individual, "_karyotype.tsv"))
cat("Karyotype inference complete:\n")
print(result)