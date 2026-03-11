#!/usr/bin/env Rscript
# =============================================================================
# infer_karyotype.R
# Infers karyotype from mosdepth 1Mb windowed coverage (CHM13v2.0)
#
# Approach:
#   - Computes chrX/autosome and chrY/autosome coverage ratios
#   - Rounds ratios to nearest 0.5 to infer integer copy counts (BeXY-inspired)
#   - Flags X chromosome abnormalities with confidence levels and reasoning
#   - Masks PAR1, PAR2, and centromere regions before arm-level comparisons
#   - Produces a coverage ratio bar plot for visual QC
#
# Literature:
#   - BeXY (Caduff et al. 2024, Genome Biol Evol): Bayesian sex karyotype inference
#     via coverage ratio rounding to nearest 0.5
#   - Ancient DNA karyotyping (Anastasiadou et al. 2024, Commun Biol):
#     mosaic Turner, XXY, XYY detected from coverage ratios
#   - SUMMER (2025, Funct Integr Genomics): clinical nanopore karyotyping pipeline
#     validates chrX/autosome ratio in males = ~0.5 of autosomes
#   - Sex chromosome aneuploidy review (Gravholt et al. 2024):
#     XXX, XYY, XXY affect ~1 in 400 live births, 90% undiagnosed
#   - CHM13v2 PAR coordinates: Nurk et al. (2022) Science 376:44-53
#   - Centromere masking: Altemose et al. (2022) Science 376:eabl9encoded
#
# Usage:
#   Rscript infer_karyotype.R <bed_gz> <summary_txt> <individual> <output_tsv>
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

args        <- commandArgs(trailingOnly = TRUE)
bed_gz      <- args[1]
summary_txt <- args[2]
individual  <- args[3]
output_tsv  <- args[4]

cat("=== infer_karyotype.R ===\n")
cat(sprintf("Individual : %s\n", individual))
cat(sprintf("BED file   : %s\n", bed_gz))
cat(sprintf("Summary    : %s\n", summary_txt))
cat("\n")

# =============================================================================
# EXCLUSION ZONES — CHM13v2.0 chrX
# PAR1/PAR2  : Nurk et al. (2022) Science 376:44-53
# Centromere : Altemose et al. (2022) Science 376:eabl9encoded
# =============================================================================
EXCL <- list(
  PAR1       = c(0,          2394410),    # maps to both X and Y → coverage doubled
  centromere = c(58000000,   62000000),   # alpha-satellite repeat array → unreliable
  PAR2       = c(153925834,  154259566)   # maps to both X and Y → coverage doubled
)

XP_START <- EXCL$PAR1[2]       + 1
XP_END   <- EXCL$centromere[1] - 1
XQ_START <- EXCL$centromere[2] + 1
XQ_END   <- EXCL$PAR2[1]       - 1

# =============================================================================
# THRESHOLDS
# =============================================================================
ARM_IMBALANCE_THRESHOLD <- 0.60   # Xp/Xq < this → arm deletion flag
MIN_AUTOSOME_COV        <- 5.0    # below this → unreliable calls

# =============================================================================
# 1. Load mosdepth summary
# =============================================================================
cat("--- Step 1: Loading mosdepth summary ---\n")

summary_df <- read_tsv(summary_txt, show_col_types = FALSE) %>%
  rename_with(tolower)

get_mean_cov <- function(chr) {
  row <- summary_df %>% filter(chrom == chr)
  if (nrow(row) == 0) {
    cat(sprintf("  WARNING: %s not found in summary\n", chr))
    return(NA_real_)
  }
  return(row$mean[1])
}

chrX_cov  <- get_mean_cov("chrX")
chrY_cov  <- get_mean_cov("chrY")
chr18_cov <- get_mean_cov("chr18")
chr21_cov <- get_mean_cov("chr21")
autosome_cov <- mean(c(chr18_cov, chr21_cov), na.rm = TRUE)

cat(sprintf("  chrX  : %.3fx\n", chrX_cov))
cat(sprintf("  chrY  : %.3fx\n", chrY_cov))
cat(sprintf("  chr18 : %.3fx\n", chr18_cov))
cat(sprintf("  chr21 : %.3fx\n", chr21_cov))
cat(sprintf("  auto  : %.3fx\n", autosome_cov))
cat("\n")

# =============================================================================
# 2. Karyotype inference — BeXY-inspired rounding approach
#    Caduff et al. (2024, Genome Biol Evol)
#    Round chrX/A and chrY/A to nearest 0.5 = integer copy count
#
#    Expected ratios:
#      chrX/A : XO=0.5, XY=0.5, XX=1.0, XXY=1.0, XXX=1.5, XXYY=1.0
#      chrY/A : no Y=0.0, Y present=0.5, YY=1.0
# =============================================================================
cat("--- Step 2: Karyotype inference (BeXY-inspired rounding) ---\n")

chrX_ratio <- chrX_cov / autosome_cov
chrY_ratio <- chrY_cov / autosome_cov

# Round to nearest 0.5 to get copy number estimates
x_copies <- round(chrX_ratio / 0.5) * 0.5
y_copies <- round(chrY_ratio / 0.5) * 0.5

cat(sprintf("  chrX/A raw   : %.3f → %.1f X copies\n", chrX_ratio, x_copies / 0.5))
cat(sprintf("  chrY/A raw   : %.3f → %.1f Y copies\n", chrY_ratio, y_copies / 0.5))

# Map copy counts to karyotype string
karyotype <- case_when(
  x_copies == 0.5 & y_copies == 0.0 ~ "XO",
  x_copies == 0.5 & y_copies == 0.5 ~ "XY",
  x_copies == 1.0 & y_copies == 0.0 ~ "XX",
  x_copies == 1.0 & y_copies == 0.5 ~ "XXY",
  x_copies == 1.5 & y_copies == 0.0 ~ "XXX",
  x_copies == 1.0 & y_copies == 1.0 ~ "XYY",
  x_copies == 1.5 & y_copies == 0.5 ~ "XXXY",
  x_copies == 1.0 & y_copies == 1.5 ~ "XYYY",
  x_copies == 1.5 & y_copies == 1.0 ~ "XXYY",
  TRUE ~ sprintf("unknown(X=%.1f,Y=%.1f)", x_copies / 0.5, y_copies / 0.5)
)

# Confidence based on autosomal coverage
karyotype_conf <- case_when(
  autosome_cov >= 20 ~ "HIGH",
  autosome_cov >= 5  ~ "MEDIUM",
  TRUE               ~ "LOW"
)

cat(sprintf("  Karyotype    : %s (%s confidence)\n", karyotype, karyotype_conf))
cat(sprintf("  Method: BeXY-inspired rounding (Caduff et al. 2024, Genome Biol Evol)\n"))
cat("\n")

# =============================================================================
# 3. Load 1Mb windowed chrX coverage and apply exclusion masks
# =============================================================================
cat("--- Step 3: Loading windowed coverage + masking ---\n")

bed <- read_tsv(bed_gz,
                col_names = c("chrom", "start", "end", "coverage"),
                show_col_types = FALSE) %>%
  filter(chrom == "chrX")

is_excluded <- function(start, end) {
  for (zone in EXCL) {
    if (start < zone[2] && end > zone[1]) return(TRUE)
  }
  return(FALSE)
}

bed <- bed %>%
  rowwise() %>%
  mutate(excluded = is_excluded(start, end)) %>%
  ungroup()

xp_windows <- bed %>% filter(!excluded, start >= XP_START, end <= XP_END)
xq_windows <- bed %>% filter(!excluded, start >= XQ_START, end <= XQ_END)

xp_mean   <- mean(xp_windows$coverage, na.rm = TRUE)
xq_mean   <- mean(xq_windows$coverage, na.rm = TRUE)
arm_ratio <- xp_mean / xq_mean

cat(sprintf("  Total chrX windows  : %d\n", nrow(bed)))
cat(sprintf("  Excluded windows    : %d (PAR1, centromere, PAR2)\n", sum(bed$excluded)))
cat(sprintf("  Xp mean (masked)    : %.3fx  (%d windows)\n", xp_mean, nrow(xp_windows)))
cat(sprintf("  Xq mean (masked)    : %.3fx  (%d windows)\n", xq_mean, nrow(xq_windows)))
cat(sprintf("  Xp/Xq ratio         : %.3f\n", arm_ratio))
cat("\n")

# =============================================================================
# 4. Generate flags with confidence + reasoning
#    NOTE: flagging logic is under review — literature review in progress
# =============================================================================
cat("--- Step 4: Flagging abnormalities ---\n")

flags    <- list()
flag_idx <- 0

add_flag <- function(name, confidence, reason, literature) {
  flag_idx <<- flag_idx + 1
  flags[[flag_idx]] <<- list(
    name       = name,
    confidence = confidence,
    reason     = reason,
    literature = literature
  )
  cat(sprintf("  [%s] %s\n    %s\n", confidence, name, reason))
}

# Flag: low coverage
if (!is.na(autosome_cov) && autosome_cov < MIN_AUTOSOME_COV) {
  add_flag(
    name       = "low_coverage_warning",
    confidence = "HIGH",
    reason     = sprintf(
      "Autosomal coverage %.2fx below 5x minimum — all calls unreliable",
      autosome_cov),
    literature = "Stancu et al. (2021) Nat Comms: minimum 5x coverage for reliable CNV calling from long reads"
  )
}

# Flag: mosaicism
if (!is.na(chrX_ratio) && chrX_ratio >= 0.60 && chrX_ratio < 0.75) {
  add_flag(
    name       = "mosaic_X_suspected",
    confidence = "MEDIUM",
    reason     = sprintf(
      "chrX/A ratio=%.3f in ambiguous zone (0.60-0.75) between XO and XX — consistent with mosaic Turner syndrome (mix of XO and XX cells)",
      chrX_ratio),
    literature = "Anastasiadou et al. (2024) Commun Biol: mosaic Turner syndrome detected from coverage ratios; Gravholt et al. (2024): 15-25% of Turner cases are mosaic 45,X/46,XX"
  )
}

# Flag: isochromosome Xq
if (!is.na(xp_mean) && !is.na(xq_mean) && xp_mean < (xq_mean * 0.25)) {
  add_flag(
    name       = "iso_Xq_suspected",
    confidence = "MEDIUM",
    reason     = sprintf(
      "Xp mean=%.3fx is less than 25%% of Xq mean=%.3fx — consistent with isochromosome Xq (Xp lost, Xq duplicated). PAR1 and centromere excluded. Directly affects X inactivation skewness.",
      xp_mean, xq_mean),
    literature = "Gravholt et al. (2024) Nat Rev Endocrinol: iso(Xq) in ~10-12% of Turner syndrome patients, detectable by Xp coverage loss and Xq gain from sequencing"
  )
}

# Flag: Xp arm deletion
if (!is.na(arm_ratio) && arm_ratio < ARM_IMBALANCE_THRESHOLD && xp_mean >= (xq_mean * 0.25)) {
  add_flag(
    name       = "Xp_deletion_suspected",
    confidence = "MEDIUM",
    reason     = sprintf(
      "Xp/Xq ratio=%.3f below threshold %.2f. Xp mean=%.3fx, Xq mean=%.3fx. PAR1 (chrX:0-2,394,410) and centromere (chrX:58-62Mb) excluded.",
      arm_ratio, ARM_IMBALANCE_THRESHOLD, xp_mean, xq_mean),
    literature = "Nurk et al. (2022) Science: PAR masking essential for arm-level coverage in CHM13v2; SUMMER (2025) Funct Integr Genomics: arm-level coverage validated in clinical nanopore pipeline"
  )
}

# Flag: Xq arm deletion
if (!is.na(arm_ratio) && arm_ratio > (1 / ARM_IMBALANCE_THRESHOLD)) {
  add_flag(
    name       = "Xq_deletion_suspected",
    confidence = "MEDIUM",
    reason     = sprintf(
      "Xp/Xq ratio=%.3f above threshold %.2f. Xp mean=%.3fx, Xq mean=%.3fx. PAR2 (chrX:153,925,834-154,259,566) and centromere (chrX:58-62Mb) excluded.",
      arm_ratio, 1 / ARM_IMBALANCE_THRESHOLD, xp_mean, xq_mean),
    literature = "Nurk et al. (2022) Science: PAR masking essential for arm-level coverage in CHM13v2; SUMMER (2025) Funct Integr Genomics: arm-level coverage validated in clinical nanopore pipeline"
  )
}

if (length(flags) == 0) cat("  No flags raised\n")
cat("\n")

# =============================================================================
# 5. Coverage ratio plot
# =============================================================================
cat("--- Step 5: Generating coverage ratio plot ---\n")

plot_data <- data.frame(
  chromosome = c("chr18", "chr21", "chrX", "chrY"),
  mean_cov   = c(chr18_cov, chr21_cov, chrX_cov, chrY_cov),
  type       = c("Autosome", "Autosome", "Sex chromosome", "Sex chromosome")
) %>%
  mutate(
    ratio      = mean_cov / autosome_cov,
    chromosome = factor(chromosome, levels = c("chr18", "chr21", "chrX", "chrY"))
  )

expected_lines <- data.frame(
  ratio = c(0.5, 1.0, 1.5),
  label = c("0.5x (1 copy)", "1.0x (2 copies)", "1.5x (3 copies)")
)

p <- ggplot(plot_data, aes(x = chromosome, y = ratio, fill = type)) +
  geom_col(width = 0.6, colour = "grey30", linewidth = 0.3) +
  geom_hline(data = expected_lines,
             aes(yintercept = ratio, linetype = label),
             colour = "grey40", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2fx", ratio)),
            vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Autosome" = "#8FBCDB", "Sex chromosome" = "#E07B5A")) +
  scale_linetype_manual(values = c("0.5x (1 copy)"   = "dotted",
                                   "1.0x (2 copies)" = "dashed",
                                   "1.5x (3 copies)" = "longdash")) +
  scale_y_continuous(limits = c(0, max(plot_data$ratio, 1.8, na.rm = TRUE) + 0.3),
                     breaks = seq(0, 2, by = 0.25)) +
  labs(
    title    = sprintf("Coverage ratio plot — %s", individual),
    subtitle = sprintf("Karyotype: %s (%s confidence) | Autosomal baseline: %.1fx",
                       karyotype, karyotype_conf, autosome_cov),
    x        = "Chromosome",
    y        = "Coverage ratio (relative to autosomal mean)",
    fill     = "Chromosome type",
    linetype = "Expected copy number"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    legend.position  = "bottom",
    legend.box       = "vertical",
    panel.grid.minor = element_blank()
  )

plot_file <- sprintf("%s_karyotype_coverage.png", individual)
ggsave(plot_file, p, width = 7, height = 5, dpi = 150)
cat(sprintf("  Plot saved: %s\n\n", plot_file))

# =============================================================================
# 6. Serialise flags
# =============================================================================
if (length(flags) > 0) {
  flag_names  <- paste(sapply(flags, `[[`, "name"),       collapse = "|")
  flag_conf   <- paste(sapply(flags, `[[`, "confidence"), collapse = "|")
  flag_reason <- paste(sapply(flags, `[[`, "reason"),     collapse = " || ")
  flag_lit    <- paste(sapply(flags, `[[`, "literature"), collapse = " || ")
} else {
  flag_names  <- "none"
  flag_conf   <- "none"
  flag_reason <- "none"
  flag_lit    <- "none"
}

# =============================================================================
# 7. Write output TSV
# =============================================================================
cat("--- Step 6: Writing output ---\n")

result <- data.frame(
  individual      = individual,
  karyotype       = karyotype,
  karyotype_conf  = karyotype_conf,
  x_copies        = x_copies / 0.5,
  y_copies        = y_copies / 0.5,
  chrX_mean_cov   = round(chrX_cov, 3),
  chrY_mean_cov   = round(chrY_cov, 3),
  auto_mean_cov   = round(autosome_cov, 3),
  chrX_ratio      = round(chrX_ratio, 3),
  chrY_ratio      = round(chrY_ratio, 3),
  xp_mean_cov     = round(xp_mean, 3),
  xq_mean_cov     = round(xq_mean, 3),
  xp_xq_ratio     = round(arm_ratio, 3),
  flags           = flag_names,
  flag_confidence = flag_conf,
  flag_reason     = flag_reason,
  flag_literature = flag_lit,
  stringsAsFactors = FALSE
)

write_tsv(result, output_tsv)

cat(sprintf("  Written: %s\n", output_tsv))
cat("\n=== DONE ===\n")
cat(sprintf("  Karyotype : %s (%s)\n", karyotype, karyotype_conf))
cat(sprintf("  Flags     : %s\n", flag_names))