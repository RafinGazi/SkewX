#!/usr/bin/env python3

# =============================================================================
# infer_karyotype_v2.py
#
# Robust karyotype inference from long-read sequencing coverage.
#
# This script improves the original infer_karyotype.R implementation by adding:
#
# - robust autosomal coverage estimation
# - MAD outlier filtering
# - X chromosome arm analysis (Xp vs Xq)
# - heterozygosity analysis from VCF
# - mosaic detection
# - RX/RY clustering validation
# - confidence scoring
#
# Inputs
# ------
# mosdepth window coverage
# mosdepth summary
# VCF file
#
# Outputs
# -------
# TSV summary
# QC plots
#
# =============================================================================

# =============================================================================
# Imports
# =============================================================================

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Logging helper
# =============================================================================

def log(message):
    print(f"[infer_karyotype] {message}")


# =============================================================================
# CHM13v2 chrX exclusion zones
# =============================================================================

PAR1 = (0, 2394410)
CENTROMERE = (58000000, 62000000)
PAR2 = (153925834, 154259566)

XP_START = PAR1[1] + 1
XP_END = CENTROMERE[0] - 1

XQ_START = CENTROMERE[1] + 1
XQ_END = PAR2[0] - 1

# =============================================================================
# Thresholds
# =============================================================================

ARM_IMBALANCE_THRESHOLD = 0.60
MIN_AUTOSOME_COV = 5.0

RX_ONE_COPY = 0.65
RX_TWO_COPY = 1.35

RY_PRESENT = 0.05
RY_TWO_COPY = 0.75

# =============================================================================
# Argument parsing
# =============================================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description="Infer karyotype from mosdepth coverage"
    )

    parser.add_argument("bed_gz", help="mosdepth window coverage file")
    parser.add_argument("summary_txt", help="mosdepth summary file")
    parser.add_argument("individual", help="sample ID")
    parser.add_argument("output_tsv", help="output TSV")
    parser.add_argument("--vcf", help="VCF file for heterozygosity analysis", required=False)

    return parser.parse_args()

# =============================================================================
# Load mosdepth summary
# =============================================================================

def load_summary(summary_file):

    df = pd.read_csv(summary_file, sep="\t")

    df.columns = [c.lower() for c in df.columns]

    cov = dict(zip(df["chrom"], df["mean"]))

    chrX_cov = cov.get("chrX", np.nan)
    chrY_cov = cov.get("chrY", np.nan)
    chr18_cov = cov.get("chr18", np.nan)
    chr21_cov = cov.get("chr21", np.nan)

    autosome_vals = [v for v in [chr18_cov, chr21_cov] if not np.isnan(v)]

    if len(autosome_vals) == 0:
        log("ERROR: No autosomal coverage available")
        autosome_cov = np.nan
    elif len(autosome_vals) == 1:
        log("WARNING: Only one autosome available for baseline")
        autosome_cov = autosome_vals[0]
    else:
        autosome_cov = np.median(autosome_vals)

    return chrX_cov, chrY_cov, chr18_cov, chr21_cov, autosome_cov

# =============================================================================
# Load window coverage
# =============================================================================

def load_windows(window_file):

    df = pd.read_csv(
        window_file,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "coverage"]
    )

    return df


# =============================================================================
# Robust coverage filtering using MAD
# =============================================================================

def filter_outliers(df, column="coverage"):

    values = df[column].values

    median = np.median(values)
    mad = np.median(np.abs(values - median))

    if mad == 0:
        return df

    mask = np.abs(values - median) < (3 * mad)

    return df[mask]

# =============================================================================
# Mask problematic regions on chrX
# =============================================================================

def mask_chrX_regions(window_df):

    x_windows = window_df[window_df["chrom"] == "chrX"].copy()

    def is_excluded(start, end):

        if start < PAR1[1] and end > PAR1[0]:
            return True

        if start < CENTROMERE[1] and end > CENTROMERE[0]:
            return True

        if start < PAR2[1] and end > PAR2[0]:
            return True

        return False

    x_windows["excluded"] = x_windows.apply(
        lambda row: is_excluded(row.start, row.end),
        axis=1
    )

    x_windows = x_windows[~x_windows["excluded"]]

    return x_windows

# =============================================================================
# Apply MAD filtering to coverage windows
# =============================================================================

def filter_chrX_windows(x_windows):

    # Remove extreme coverage windows using MAD filtering
    # Helps stabilise arm-level coverage estimates

    x_windows = filter_outliers(x_windows, "coverage")

    return x_windows

# =============================================================================
# Separate X chromosome arms
# =============================================================================

def split_chrX_arms(x_windows):

    xp = x_windows[
        (x_windows["start"] >= XP_START) &
        (x_windows["end"] <= XP_END)
    ]

    xq = x_windows[
        (x_windows["start"] >= XQ_START) &
        (x_windows["end"] <= XQ_END)
    ]

    xp_mean = np.median(xp["coverage"]) if len(xp) > 0 else np.nan
    xq_mean = np.median(xq["coverage"]) if len(xq) > 0 else np.nan

    if np.isnan(xp_mean) or np.isnan(xq_mean) or xq_mean <= 0:
        arm_ratio = np.nan
    else:
        arm_ratio = xp_mean / xq_mean

    return xp_mean, xq_mean, arm_ratio

# =============================================================================
# Mosaic detection using coverage variance
# =============================================================================

def detect_mosaicism(x_windows, autosome_windows=None):

    if len(x_windows) == 0:
        return np.nan, np.nan, False, np.nan

    cov = x_windows["coverage"].values

    median_cov = np.median(cov)
    mad = np.median(np.abs(cov - median_cov))

    mad_ratio = mad / median_cov if median_cov > 0 else np.nan

    # Optional normalization using autosomes (better science)
    if autosome_windows is not None and len(autosome_windows) > 0:
        auto_cov = autosome_windows["coverage"].values
        auto_median = np.median(auto_cov)
        auto_mad = np.median(np.abs(auto_cov - auto_median))

        if auto_median > 0:
            auto_mad_ratio = auto_mad / auto_median
            z_score = mad_ratio / auto_mad_ratio if auto_mad_ratio > 0 else np.nan
        else:
            z_score = np.nan
    else:
        z_score = np.nan

    # Research-style threshold
    mosaic_flag = z_score > 2.5 if not np.isnan(z_score) else mad_ratio > 0.20

    return mad, mad_ratio, mosaic_flag, z_score

# =============================================================================
# Detect X chromosome structural abnormalities
# =============================================================================

def detect_arm_abnormalities(xp_mean, xq_mean, arm_ratio):

    flags = []

    # Case 1: Xp completely missing
    if np.isnan(xp_mean) and not np.isnan(xq_mean):
        flags.append("Xp_deletion")
        return flags

    # Case 2: Xq completely missing
    if np.isnan(xq_mean) and not np.isnan(xp_mean):
        flags.append("Xq_deletion")
        return flags

    # If both missing → cannot infer
    if np.isnan(xp_mean) and np.isnan(xq_mean):
        return flags

    # Normal imbalance detection
    if not np.isnan(arm_ratio):

        # Isochromosome Xq (very strong imbalance)
        if xp_mean < (0.25 * xq_mean):
            flags.append("iso_Xq")

        elif arm_ratio < ARM_IMBALANCE_THRESHOLD:
            flags.append("Xp_deletion")

        elif arm_ratio > (1 / ARM_IMBALANCE_THRESHOLD):
            flags.append("Xq_deletion")

    return flags

# =============================================================================
# Heterozygosity analysis from VCF
# =============================================================================

def compute_heterozygosity(vcf_file):

    if vcf_file is None:
        log("WARNING: No VCF provided — heterozygosity skipped")
        return np.nan, 0, 0

    het_count = 0
    total_count = 0

    with open(vcf_file, "r") as f:
        for line in f:

            # Skip headers
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")

            chrom = fields[0]

            # Only chrX
            if chrom not in ["chrX", "X"]:
                continue
            
            if len(fields) < 10:
                continue

            format_fields = fields[8].split(":")
            sample_fields = fields[9].split(":")

            if len(format_fields) != len(sample_fields):
                continue

            format_dict = dict(zip(format_fields, sample_fields))

            genotype = format_dict.get("GT", "./.")
            depth_val = format_dict.get("DP", "0")

            try:
                depth = int(depth_val)
            except:
                depth = 0

            # Skip missing genotype
            if genotype in ["./.", ".|."]:
                continue

            # Skip low-quality variants
            if depth < 10:
                continue

            # Accept both unphased and phased genotypes
            valid_genotypes = ["0/0", "0/1", "1/0", "1/1", "0|1", "1|0"]
            het_genotypes = {"0/1","1/0","0|1","1|0"}

            if genotype in valid_genotypes:
                total_count += 1

                if genotype in het_genotypes:
                    het_count += 1  

    hx = het_count / total_count if total_count > 0 else np.nan

    return hx, het_count, total_count

def interpret_heterozygosity(hx):

    if np.isnan(hx):
        return "unknown"

    if hx < 0.05:
        return "low"
    elif hx < 0.25:
        return "moderate"
    else:
        return "high"

# =============================================================================
# Coverage ratio calculations
# =============================================================================

def compute_ratios(chrX_cov, chrY_cov, autosome_cov):

    if autosome_cov <= 0 or np.isnan(autosome_cov):
        log("ERROR: Invalid autosome coverage")
        return np.nan, np.nan

    # Handle missing chrY (VERY IMPORTANT)
    if np.isnan(chrY_cov):
        log("WARNING: chrY coverage missing — assuming 0")
        chrY_cov = 0.0

    # chrX missing is still a real error
    if np.isnan(chrX_cov):
        log("ERROR: Missing chrX coverage")
        return np.nan, np.nan   

    rx = chrX_cov / autosome_cov
    # Normalize RX using diploid expectation (~XX samples)
    rx_normalized = rx / 2.0
    ry = chrY_cov / autosome_cov
    rx = rx_normalized

    if np.isnan(rx) or np.isnan(ry):
        log("ERROR: Invalid RX/RY values — check input coverage")

    return rx, ry

# =============================================================================
# Infer chromosome copy numbers
# =============================================================================

def infer_copy_numbers(rx, ry):

    # Infer X copies
    if rx < RX_ONE_COPY:
        x_copies = 1
    elif rx < RX_TWO_COPY:
        x_copies = 2
    else:
        x_copies = 3

    # Infer Y copies
    if ry < RY_PRESENT:
        y_copies = 0
    elif ry < RY_TWO_COPY:
        y_copies = 1
    else:
        y_copies = 2

    return x_copies, y_copies


# =============================================================================
# Convert copy numbers to karyotype label
# =============================================================================

def determine_karyotype(x_copies, y_copies):

    mapping = {
        (2,0): "XX",
        (1,1): "XY",
        (1,0): "XO",
        (2,1): "XXY",
        (3,0): "XXX",
        (1,2): "XYY",
        (3,1): "XXXY"
    }

    return mapping.get((x_copies, y_copies), "unknown")

# =============================================================================
# Coverage ratio plot
# =============================================================================

def plot_coverage_ratios(chr18_cov, chr21_cov, chrX_cov, chrY_cov, autosome_cov, sample):

    if autosome_cov <= 0 or np.isnan(autosome_cov):
        log("WARNING: Skipping coverage ratio plot due to invalid autosome coverage")
        return

    ratios = [
        chr18_cov / autosome_cov,
        chr21_cov / autosome_cov,
        chrX_cov / autosome_cov,
        chrY_cov / autosome_cov
    ]

    labels = ["chr18", "chr21", "chrX", "chrY"]

    # Autosomes = blue, sex chromosomes = orange
    colors = ["#8FBCDB", "#8FBCDB", "#E07B5A", "#E07B5A"]

    plt.figure(figsize=(7,5))

    plt.bar(labels, ratios, color=colors, edgecolor="black")

    # Expected copy-number reference lines
    plt.axhline(0.5, linestyle="dotted", color="grey", label="1 copy")
    plt.axhline(1.0, linestyle="dashed", color="grey", label="2 copies")
    plt.axhline(1.5, linestyle="dotted", color="grey", label="3 copies")

    plt.ylabel("Coverage ratio (relative to autosomes)")
    plt.title(f"Coverage ratios — {sample}")

    valid_ratios = [r for r in ratios if not np.isnan(r)]

    if len(valid_ratios) == 0:
        log("WARNING: No valid ratios for plotting")
        return
    plt.ylim(0, max(valid_ratios) + 0.4)

    plt.legend()    

    plt.savefig(f"{sample}_coverage_ratios.png", dpi=150)
    plt.close()

# =============================================================================
# RX vs RY plot
# =============================================================================

def plot_rx_ry(rx, ry, sample):

    plt.figure(figsize=(5,5))

    plt.scatter(rx, ry, color="red", s=80)

    plt.xlabel("RX (X coverage / autosome coverage)")
    plt.ylabel("RY (Y coverage / autosome coverage)")

    plt.title(f"RX vs RY — {sample}")

    plt.axvline(0.5, linestyle="dotted", color="grey")
    plt.axvline(1.0, linestyle="dotted", color="grey")

    plt.axhline(0.5, linestyle="dotted", color="grey")

    plt.xlim(0, 2)
    plt.ylim(0, 2)

    plt.text(1.0, 0.05, "XX", fontsize=10)
    plt.text(0.5, 0.5, "XY", fontsize=10)
    plt.text(0.5, 0.05, "XO", fontsize=10)

    plt.savefig(f"{sample}_rx_ry.png", dpi=150)
    plt.close()

# =============================================================================
# chrX coverage profile plot
# =============================================================================

def plot_chrX_coverage(x_windows, sample):

    if len(x_windows) == 0:
        log("WARNING: No chrX windows — skipping coverage plot")
        return

    plt.figure(figsize=(10,4))

    plt.scatter(
        x_windows["start"],
        x_windows["coverage"],
        s=10,
        alpha=0.7
    )

    plt.xlabel("Genomic position on chrX")
    plt.ylabel("Coverage")

    plt.title(f"chrX coverage profile — {sample}")

    plt.savefig(f"{sample}_chrX_coverage.png", dpi=150)
    plt.close()

# =============================================================================
# Confidence scoring
# =============================================================================
def compute_confidence(rx, ry, hx, mosaic_flag):

    if np.isnan(rx) or np.isnan(ry):
        return "LOW"

    score = 0

    # Coverage confidence
    if not np.isnan(rx) and not np.isnan(ry):
        score += 1

    # Heterozygosity support
    if not np.isnan(hx):
        score += 1

    # Penalize mosaic
    if mosaic_flag:
        score -= 1

    if score >= 2:
        return "HIGH"
    elif score == 1:
        return "MEDIUM"
    else:
        return "LOW"

# =============================================================================
# Cross-check consistency between signals
# =============================================================================

def check_consistency(karyotype, hx, mosaic_flag):

    flags = []

    # Expected heterozygosity patterns
    if karyotype == "XX" and not np.isnan(hx):
        if hx < 0.08:
            flags.append("low_het_for_XX")

    elif karyotype in ["XY", "XO"]:
        if hx > 0.1:
            flags.append("unexpected_heterozygosity")

    elif karyotype in ["XXY", "XXX"]:
        if hx < 0.1:
            flags.append("low_het_for_extra_X")

    # Mosaic override
    if mosaic_flag:
        flags.append("mosaic_sample")

    status = "consistent" if len(flags) == 0 else "inconsistent"

    return status, flags

# =============================================================================
# Main execution
# =============================================================================

def main():

    args = parse_args()

    log(f"Sample: {args.individual}")

    # Load mosdepth summary
    log("Loading mosdepth summary")
    chrX_cov, chrY_cov, chr18_cov, chr21_cov, autosome_cov = load_summary(args.summary_txt)

    log(f"chrX coverage: {chrX_cov}")
    log(f"chrY coverage: {chrY_cov}")
    log(f"autosome baseline: {autosome_cov}")

    if np.isnan(autosome_cov):
        log("ERROR: Autosome coverage is NaN — writing fallback output")

        result = pd.DataFrame([{
            "individual": args.individual,
            "karyotype": "unknown",
            "karyotype_confidence": "LOW",
            "chrX_ratio": np.nan,
            "chrY_ratio": np.nan,
            "cohort_qc_flag": "pending",
            "consistency_status": "unknown",
            "consistency_flags": "autosome_missing"
        }])

        result.to_csv(args.output_tsv, sep="\t", index=False)
        return

    elif autosome_cov < MIN_AUTOSOME_COV:
        log("WARNING: Low autosomal coverage — results may be unreliable")

    # Load window coverage
    log("Loading window coverage")
    windows = load_windows(args.bed_gz)

    # Process chrX windows
    log("Masking PAR and centromere regions")
    x_windows = mask_chrX_regions(windows)

    log("Applying MAD filtering")
    x_windows = filter_chrX_windows(x_windows)
    if len(x_windows) < 50:
        log("WARNING: Very few chrX windows — results may be unstable")

    log("Generating chrX coverage profile plot")
    plot_chrX_coverage(x_windows, args.individual)

    log("Computing Xp/Xq arm coverage")
    xp_mean, xq_mean, arm_ratio = split_chrX_arms(x_windows)

    log("Detecting mosaicism from coverage variance")

    autosome_windows = windows[windows["chrom"].isin(["chr18", "chr21"])]
    if len(autosome_windows) == 0:
        log("WARNING: No autosomal windows available for mosaic normalization")
    autosome_windows = filter_outliers(autosome_windows, "coverage")
    mad, mad_ratio, mosaic_flag, z_score = detect_mosaicism(x_windows, autosome_windows)

    log(f"chrX MAD: {mad}")
    log(f"MAD/coverage ratio: {mad_ratio}")
    log(f"Mosaic flag: {mosaic_flag}")
    log(f"Mosaic z-score: {z_score}")

    log(f"Xp mean: {xp_mean}")
    log(f"Xq mean: {xq_mean}")
    log(f"Xp/Xq ratio: {arm_ratio}")

    mosaic_status = "mosaic_candidate" if mosaic_flag else "no_mosaic_signal"

    log("Detecting arm-level abnormalities")

    arm_flags = detect_arm_abnormalities(xp_mean, xq_mean, arm_ratio)

    log(f"Arm-level flags: {arm_flags if arm_flags else 'none'}")

    # Coverage ratios
    rx, ry = compute_ratios(chrX_cov, chrY_cov, autosome_cov)

    log(f"RX: {rx}")
    log(f"RY: {ry}")

    if not np.isnan(rx) and not np.isnan(ry):
        if rx > 2 or ry > 2:
            log("WARNING: Extreme RX/RY values — possible coverage issue")

    log(f"[QC] RX/RY ready for cohort comparison")

    # Infer copy numbers
    if np.isnan(rx) or np.isnan(ry):
        log("ERROR: Cannot infer karyotype due to invalid ratios")
        karyotype = "unknown"
        x_copies, y_copies = np.nan, np.nan
    else:
        x_copies, y_copies = infer_copy_numbers(rx, ry)
        # Determine karyotype
        karyotype = determine_karyotype(x_copies, y_copies)

    log(f"Inferred X copies: {x_copies}")
    log(f"Inferred Y copies: {y_copies}")
    log(f"Inferred karyotype: {karyotype}")

    log("Computing chrX heterozygosity")

    hx, het_count, total_count = compute_heterozygosity(args.vcf)
    hx_label = interpret_heterozygosity(hx)

    log(f"HX (heterozygosity): {hx}")
    log(f"Heterozygous variants: {het_count}")
    log(f"Total variants: {total_count}")

    confidence = compute_confidence(rx, ry, hx, mosaic_flag)
    log(f"Karyotype confidence: {confidence}")

    log("Performing cross-signal consistency check")

    consistency_status, consistency_flags = check_consistency(
        karyotype,
        hx,
        mosaic_flag
    )

    log(f"Consistency: {consistency_status}")
    log(f"Consistency flags: {consistency_flags if consistency_flags else 'none'}")

    log("Generating coverage ratio plot")

    plot_coverage_ratios(
        chr18_cov,
        chr21_cov,
        chrX_cov,
        chrY_cov,
        autosome_cov,
        args.individual
    )

    if not np.isnan(rx) and not np.isnan(ry):
        plot_rx_ry(rx, ry, args.individual)
    else:
        log("Skipping RX/RY plot due to invalid values")

    # Save result
    result = pd.DataFrame([{
    "individual": args.individual,
    "karyotype": karyotype,
    "karyotype_conf": confidence,
    "x_copies": int(x_copies) if not np.isnan(x_copies) else np.nan,
    "y_copies": int(y_copies) if not np.isnan(y_copies) else np.nan,
    "auto_mean_cov": round(autosome_cov, 3) if not np.isnan(autosome_cov) else np.nan,
    "flags": "none",
    "flag_confidence": "none",
    "flag_reason": "none",
    "flag_literature": "none",
    "n_chrX_windows": len(x_windows),
    "chrX_ratio": round(rx, 3) if not np.isnan(rx) else np.nan,
    "chrY_ratio": round(ry, 3) if not np.isnan(ry) else np.nan,
    "chrX_heterozygosity": round(hx, 3) if not np.isnan(hx) else np.nan,
    "chrX_het_level": hx_label,
    "chrX_het_sites": het_count,
    "chrX_total_sites": total_count,
    "consistency_status": consistency_status,
    "consistency_flags": "|".join(consistency_flags) if consistency_flags else "none",
    "xp_xq_ratio": round(arm_ratio, 3) if not np.isnan(arm_ratio) else np.nan,
    "xp_mean_cov": round(xp_mean, 3) if not np.isnan(xp_mean) else np.nan,
    "xq_mean_cov": round(xq_mean, 3) if not np.isnan(xq_mean) else np.nan,
    "chrX_mad": round(mad, 3) if not np.isnan(mad) else np.nan,
    "chrX_mad_ratio": round(mad_ratio, 3) if not np.isnan(mad_ratio) else np.nan,
    "chrX_mosaic_zscore": round(z_score, 3) if not np.isnan(z_score) else np.nan,
    "mosaic_status": mosaic_status,
    "arm_flags": "|".join(arm_flags) if arm_flags else "none",
    "cohort_qc_flag": "pending",
}])

    log("Writing output TSV")
    result.to_csv(args.output_tsv, sep="\t", index=False)

    log("Karyotype inference complete — awaiting cohort QC")
    log("Done")

    


if __name__ == "__main__":
    main()

