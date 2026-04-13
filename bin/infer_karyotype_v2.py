#!/usr/bin/env python3

# =============================================================================
# infer_karyotype.py
#
# Simple karyotype inference from mosdepth coverage.
#
# Logic:
#   1. Read mosdepth summary → compute RX, RY → infer karyotype
#   2. Load windows for all samples (needed for genome coverage plot)
#   3. Run arm analysis for all karyotypes except XY and XYY
#   3b. Rescue XO → XX if arm flags detected (structural deletion, not monosomy)
#   4. If not XX → write TSV + plots, exit code 2
#   5. If XX → write TSV + all plots, exit code 0
#
# Inputs:
#   mosdepth regions bed.gz   (window coverage)
#   mosdepth summary.txt      (per-chromosome mean coverage)
#   sample ID
#   output TSV path
#
# Outputs:
#   TSV with karyotype call and flags
#   PNG plots
# =============================================================================

import argparse
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Thresholds
# =============================================================================

RX_ONE_COPY  = 0.65
RX_TWO_COPY  = 1.35
RY_PRESENT   = 0.05
RY_TWO_COPY  = 0.75

ARM_DELETION_THRESHOLD = 0.60
ARM_PARTIAL_THRESHOLD  = 0.75

# Karyotypes where arm analysis is not meaningful
SKIP_ARM = {"XY", "XYY"}

# CHM13v2 chrX arm boundaries (split Xp vs Xq only, no masking)
XP_START =   2_394_411
XP_END   =  57_999_999
XQ_START =  62_000_001
XQ_END   = 153_925_833

# =============================================================================
# Argument parsing
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Infer karyotype from mosdepth coverage"
    )
    parser.add_argument("bed_gz",      help="mosdepth window coverage (.regions.bed.gz)")
    parser.add_argument("summary_txt", help="mosdepth summary file")
    parser.add_argument("individual",  help="sample ID")
    parser.add_argument("output_tsv",  help="output TSV path")
    return parser.parse_args()

# =============================================================================
# Load mosdepth summary
# =============================================================================

def load_summary(summary_file):
    df = pd.read_csv(summary_file, sep="\t")
    df.columns = [c.lower() for c in df.columns]
    cov = dict(zip(df["chrom"], df["mean"]))

    chrX_cov = cov.get("chrX", np.nan)
    chrY_cov = cov.get("chrY", 0.0)

    autosome_vals = [v for k, v in cov.items() if k in ("chr18", "chr21")]
    autosome_cov  = np.nanmean(autosome_vals) if autosome_vals else np.nan

    chr18_cov = cov.get("chr18", np.nan)
    chr21_cov = cov.get("chr21", np.nan)

    return chrX_cov, chrY_cov, chr18_cov, chr21_cov, autosome_cov

# =============================================================================
# Compute RX / RY
# =============================================================================

def compute_ratios(chrX_cov, chrY_cov, autosome_cov):
    rx = chrX_cov / autosome_cov
    ry = chrY_cov / autosome_cov
    return rx, ry

# =============================================================================
# Infer karyotype from RX / RY
# =============================================================================

def infer_karyotype(rx, ry):
    x = 1 if rx < RX_ONE_COPY else (2 if rx < RX_TWO_COPY else 3)
    y = 0 if ry < RY_PRESENT  else (1 if ry < RY_TWO_COPY  else 2)

    mapping = {
        (2, 0): "XX",
        (1, 1): "XY",
        (1, 0): "XO",
        (2, 1): "XXY",
        (3, 0): "XXX",
        (1, 2): "XYY",
        (3, 1): "XXXY",
    }
    return mapping.get((x, y), "unknown")

# =============================================================================
# Load window coverage
# =============================================================================

def load_windows(bed_gz):
    df = pd.read_csv(
        bed_gz, sep="\t", header=None,
        names=["chrom", "start", "end", "coverage"]
    )
    return df

# =============================================================================
# Xp / Xq arm analysis
#
# Uses absolute coverage thresholds (relative to autosome_cov) to detect
# full arm deletions where one arm is near zero — these break ratio-based
# detection because the denominator collapses. Ratio-based logic is kept
# as a fallback for partial deletions.
# =============================================================================

def arm_analysis(windows, autosome_cov):
    x = windows[windows["chrom"] == "chrX"].copy()

    xp = x[(x["start"] >= XP_START) & (x["end"] <= XP_END)]
    xq = x[(x["start"] >= XQ_START) & (x["end"] <= XQ_END)]

    xp_mean = np.median(xp["coverage"]) if len(xp) > 0 else np.nan
    xq_mean = np.median(xq["coverage"]) if len(xq) > 0 else np.nan

    # Normalise relative to autosome
    xp_norm = xp_mean / autosome_cov if autosome_cov > 0 else np.nan
    xq_norm = xq_mean / autosome_cov if autosome_cov > 0 else np.nan

    # Ratio (for reporting only)
    if np.isnan(xp_mean) or np.isnan(xq_mean) or xq_mean == 0:
        arm_ratio = np.nan
    else:
        arm_ratio = xp_mean / xq_mean

    flags = []

    # --- Key idea: detect imbalance, not absolute drop ---
    if not np.isnan(xp_norm) and not np.isnan(xq_norm):

        arm_diff = abs(xp_norm - xq_norm)

        DIFF_THRESHOLD = 0.25     # detects imbalance
        ZERO_THRESHOLD = 0.2      # near-zero → full deletion
        PARTIAL_THRESHOLD = 0.75  # reduced but not zero

        if arm_diff > DIFF_THRESHOLD:

            # Xp affected
            if xp_norm < xq_norm:
                if xp_norm < ZERO_THRESHOLD:
                    flags.append("Xp_deletion")
                elif xp_norm < PARTIAL_THRESHOLD:
                    flags.append("Xp_partial_deletion")

            # Xq affected
            elif xq_norm < xp_norm:
                if xq_norm < ZERO_THRESHOLD:
                    flags.append("Xq_deletion")
                elif xq_norm < PARTIAL_THRESHOLD:
                    flags.append("Xq_partial_deletion")

    return xp_mean, xq_mean, arm_ratio, flags

# =============================================================================
# Plots
# =============================================================================

def plot_coverage_ratios(chr18_cov, chr21_cov, chrX_cov, chrY_cov,
                         autosome_cov, sample):
    labels = ["chr18", "chr21", "chrX", "chrY"]
    ratios = [
        chr18_cov / autosome_cov,
        chr21_cov / autosome_cov,
        chrX_cov  / autosome_cov,
        chrY_cov  / autosome_cov,
    ]
    colors = ["#8FBCDB", "#8FBCDB", "#E07B5A", "#E07B5A"]

    plt.figure(figsize=(7, 5))
    plt.bar(labels, ratios, color=colors, edgecolor="black")
    plt.axhline(0.5, linestyle="dotted", color="grey", label="1 copy (0.5)")
    plt.axhline(1.0, linestyle="dashed", color="grey", label="2 copies (1.0)")
    plt.axhline(1.5, linestyle="dotted", color="grey", label="3 copies (1.5)")
    plt.ylabel("Coverage ratio (relative to autosome mean)")
    plt.title(f"Coverage ratios — {sample}")
    plt.ylim(0, max([r for r in ratios if not np.isnan(r)]) + 0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{sample}_coverage_ratios.png", dpi=150)
    plt.close()


def plot_genome_coverage(windows, autosome_cov, sample):
    chroms = ["chr18", "chr21", "chrX", "chrY"]
    colors = ["#8FBCDB", "#8FBCDB", "#E07B5A", "#E07B5A"]

    all_cov = windows[windows["chrom"].isin(chroms)]["coverage"]
    y_max = all_cov.max() * 1.15 if len(all_cov) > 0 else autosome_cov * 2

    fig, axes = plt.subplots(1, 4, figsize=(18, 4), sharey=True)
    fig.suptitle(f"Genome coverage profile — {sample}", fontsize=13)

    for ax, chrom, color in zip(axes, chroms, colors):
        sub = windows[windows["chrom"] == chrom]

        if len(sub) > 0:
            ax.scatter(
                sub["start"] / 1e6,
                sub["coverage"],
                s=15, alpha=0.7, color=color
            )
        else:
            ax.text(0.5, 0.5, "no data", ha="center", va="center",
                    transform=ax.transAxes, color="grey")

        ax.axhline(autosome_cov, linestyle="dashed", color="orange",
                   linewidth=1.2, label=f"Autosome mean ({autosome_cov:.1f}x)")
        ax.axhline(autosome_cov / 2, linestyle="dotted", color="grey",
                   linewidth=1.0, label=f"Half autosome ({autosome_cov/2:.1f}x)")

        ax.set_title(chrom)
        ax.set_xlabel("Position (Mb)")
        ax.set_xlim(0, None)
        ax.set_ylim(0, y_max)

    axes[0].set_ylabel("Coverage")
    axes[0].legend(fontsize=7)

    plt.tight_layout()
    plt.savefig(f"{sample}_genome_coverage.png", dpi=150)
    plt.close()


def plot_xp_xq(xp_mean, xq_mean, autosome_cov, sample):
    if np.isnan(xp_mean) or np.isnan(xq_mean):
        return

    plt.figure(figsize=(5, 5))
    plt.bar(["Xp", "Xq"], [xp_mean, xq_mean],
            color=["#6BAED6", "#2171B5"], edgecolor="black")
    plt.axhline(autosome_cov, linestyle="dashed", color="orange",
                linewidth=1.5, label=f"Autosome mean ({autosome_cov:.1f}x)")
    plt.axhline(autosome_cov / 2, linestyle="dotted", color="grey",
                linewidth=1.2, label=f"Half autosome ({autosome_cov/2:.1f}x)")
    plt.ylabel("Median coverage")
    plt.title(f"Xp vs Xq arm coverage — {sample}")
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(f"{sample}_xp_xq.png", dpi=150)
    plt.close()


def plot_rx_ry(rx, ry, sample):
    plt.figure(figsize=(5, 5))
    plt.scatter(rx, ry, color="red", s=100, zorder=5)
    plt.axvline(0.5, linestyle="dotted", color="grey")
    plt.axvline(1.0, linestyle="dotted", color="grey")
    plt.axhline(0.5, linestyle="dotted", color="grey")
    for (label, x_pos, y_pos) in [
        ("XX",  1.0, 0.02),
        ("XY",  0.5, 0.5),
        ("XO",  0.5, 0.02),
        ("XXY", 1.0, 0.5),
        ("XXX", 1.5, 0.02),
        ("XYY", 0.5, 1.0),
    ]:
        plt.text(x_pos, y_pos, label, fontsize=9, color="dimgrey")
    plt.xlabel("RX (chrX / autosome)")
    plt.ylabel("RY (chrY / autosome)")
    plt.title(f"RX vs RY — {sample}")
    plt.xlim(0, 2)
    plt.ylim(0, 2)
    plt.tight_layout()
    plt.savefig(f"{sample}_rx_ry.png", dpi=150)
    plt.close()

# =============================================================================
# Write output TSV
# =============================================================================

def write_output(path, individual, raw_karyotype, karyotype, status, rx, ry,
                 xp_mean, xq_mean, arm_ratio, arm_flags,
                 autosome_cov, qc_flag):
    row = {
        "individual":      individual,
        "raw_karyotype":   raw_karyotype,   # original call before rescue
        "karyotype":       karyotype,        # final call after rescue
        "status":          status,
        "chrX_ratio":      round(rx, 3)        if not np.isnan(rx)        else np.nan,
        "chrY_ratio":      round(ry, 3)        if not np.isnan(ry)        else np.nan,
        "auto_mean_cov":   round(autosome_cov, 3),
        "xp_mean_cov":     round(xp_mean, 3)   if not np.isnan(xp_mean)   else np.nan,
        "xq_mean_cov":     round(xq_mean, 3)   if not np.isnan(xq_mean)   else np.nan,
        "xp_xq_ratio":     round(arm_ratio, 3) if not np.isnan(arm_ratio) else np.nan,
        "arm_flags":       "|".join(arm_flags) if arm_flags else "none",
        "qc_flag":         qc_flag,
    }
    pd.DataFrame([row]).to_csv(path, sep="\t", index=False)

# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    sample = args.individual

    print(f"[infer_karyotype] Sample: {sample}")

    # STEP 1: Load summary and compute ratios
    chrX_cov, chrY_cov, chr18_cov, chr21_cov, autosome_cov = load_summary(args.summary_txt)

    if np.isnan(autosome_cov) or autosome_cov == 0:
        print(f"[infer_karyotype] ERROR: Cannot compute autosome baseline — check summary file")
        sys.exit(1)

    print(f"[infer_karyotype] chrX={chrX_cov:.2f}x  chrY={chrY_cov:.2f}x  autosome={autosome_cov:.2f}x")

    rx, ry = compute_ratios(chrX_cov, chrY_cov, autosome_cov)
    raw_karyotype = infer_karyotype(rx, ry)
    karyotype = raw_karyotype  # may be updated by rescue in STEP 3b

    print(f"[infer_karyotype] RX={rx:.3f}  RY={ry:.3f}  → {karyotype}")

    # STEP 2: Load windows for all samples
    print(f"[infer_karyotype] Loading windows")
    windows = load_windows(args.bed_gz)

    # STEP 3: Arm analysis — skip only for XY and XYY
    if karyotype not in SKIP_ARM:
        print(f"[infer_karyotype] Running Xp/Xq arm analysis")
        xp_mean, xq_mean, arm_ratio, arm_flags = arm_analysis(windows, autosome_cov)
        arm_ratio_str = f"{arm_ratio:.3f}" if not np.isnan(arm_ratio) else "nan"
        print(f"[infer_karyotype] Xp={xp_mean:.2f}x  Xq={xq_mean:.2f}x  ratio={arm_ratio_str}")
        if arm_flags:
            print(f"[infer_karyotype] Arm flags: {arm_flags}")
    else:
        print(f"[infer_karyotype] Skipping arm analysis — not meaningful for {karyotype}")
        xp_mean, xq_mean, arm_ratio, arm_flags = np.nan, np.nan, np.nan, []

    # STEP 3b: Structural handling
    if arm_flags:
        print(f"[infer_karyotype] Structural abnormality detected → marking as unknown")
        karyotype = "unknown"

    # STEP 4: Determine QC flag

    if arm_flags:
        qc_flag = f"flagged:{','.join(arm_flags)}"
    elif karyotype in ["XX", "XXX", "XXY"]:
        qc_flag = "pass"
    else:
        qc_flag = f"skipped:{karyotype}"

    print(f"[infer_karyotype] QC flag: {qc_flag}")

    # Derive high-level status
    if karyotype == "XX":
        status = "diploid-pass"
    elif karyotype in ["XXX", "XXY"]:
        status = "non-diploid-pass"
    elif arm_flags:
        status = "structural-flagged"
    else:
        status = "skipped"

    # STEP 5: Write TSV
    write_output(
        args.output_tsv, sample, raw_karyotype, karyotype, status, rx, ry,
        xp_mean, xq_mean, arm_ratio, arm_flags,
        autosome_cov, qc_flag
    )

    # STEP 6: Plots — all samples get these
    plot_coverage_ratios(chr18_cov, chr21_cov, chrX_cov, chrY_cov, autosome_cov, sample)
    plot_genome_coverage(windows, autosome_cov, sample)
    plot_rx_ry(rx, ry, sample)

    # Xp/Xq plot only if arm analysis was run
    if karyotype not in SKIP_ARM:
        plot_xp_xq(xp_mean, xq_mean, autosome_cov, sample)

    print(f"[infer_karyotype] Final QC: {sample} → {qc_flag}")

    print(f"[infer_karyotype] Done → {args.output_tsv}")


if __name__ == "__main__":
    main()