#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler


def log(msg):
    print(f"[cohort_qc] {msg}")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_tsvs", nargs="+")
    parser.add_argument("--out_prefix", default="cohort")
    return parser.parse_args()


# =============================================================================
# Load data
# =============================================================================
def load_data(files):
    dfs = [pd.read_csv(f, sep="\t") for f in files]
    df = pd.concat(dfs, ignore_index=True)
    df = df.dropna(subset=["chrX_ratio", "chrY_ratio"])

    # Safeguard: ensure arm columns exist if any input TSV is missing them
    for col in ["arm_flags", "arm_confidence"]:
        if col not in df.columns:
            df[col] = "MISSING"
        else:
            df[col] = df[col].fillna("MISSING")  # handle partial missing

    return df


# =============================================================================
# Expected karyotype centers
# =============================================================================
EXPECTED = {
    "XX":  (1.0, 0.0),
    "XY":  (0.5, 0.5),
    "XXX": (1.5, 0.0),
    "XXY": (1.0, 0.5),
    "XYY": (0.5, 1.0),
    "XO":  (0.5, 0.0)
}


# =============================================================================
# Distance to expected
# =============================================================================
def expected_distance(row):

    k = row["karyotype"]
    rx = row["chrX_ratio"]
    ry = row["chrY_ratio"]

    if k in EXPECTED:
        ex, ey = EXPECTED[k]
        return np.sqrt((rx - ex)**2 + (ry - ey)**2)

    return np.nan


# =============================================================================
# DBSCAN clustering
# =============================================================================
def run_clustering(df):

    X = df[["chrX_ratio", "chrY_ratio"]].values

    X_scaled = StandardScaler().fit_transform(X)

    clustering = DBSCAN(eps=0.8, min_samples=2).fit(X_scaled)

    df["cluster_id"] = clustering.labels_

    return df


# =============================================================================
# Compute cluster centers
# =============================================================================
def compute_cluster_centers(df):

    centers = {}

    for cid in df["cluster_id"].unique():

        if cid == -1:
            continue

        sub = df[df["cluster_id"] == cid]

        cx = sub["chrX_ratio"].median()
        cy = sub["chrY_ratio"].median()

        centers[cid] = (cx, cy)

    return centers


# =============================================================================
# Assign QC flags
# =============================================================================
def assign_qc(df, centers):

    qc_flags = []
    confidences = []
    mosaic_flags = []

    for _, row in df.iterrows():

        cid = row["cluster_id"]
        rx = row["chrX_ratio"]
        ry = row["chrY_ratio"]

        # distance to expected
        dist_exp = expected_distance(row)

        # noise cluster
        if cid == -1:
            qc_flags.append("OUTLIER")
            confidences.append("LOW")
            mosaic_flags.append("possible_mosaic")
            continue

        # distance to cluster center
        cx, cy = centers[cid]
        dist_cluster = np.sqrt((rx - cx)**2 + (ry - cy)**2)

        # consistency check
        consistent = dist_exp < 0.25 if not np.isnan(dist_exp) else False

        # mosaic detection
        mosaic = dist_cluster > 0.25 and dist_exp > 0.25

        # assign QC
        if not consistent:
            qc_flags.append("INCONSISTENT")
            confidences.append("LOW")

        elif mosaic:
            qc_flags.append("MOSAIC_CANDIDATE")
            confidences.append("MEDIUM")

        elif dist_cluster < 0.15:
            qc_flags.append("PASS")
            confidences.append("HIGH")

        else:
            qc_flags.append("PASS")
            confidences.append("MEDIUM")

        mosaic_flags.append("yes" if mosaic else "no")

    df["cohort_qc_flag"] = qc_flags
    df["cohort_confidence"] = confidences
    df["mosaic_flag"] = mosaic_flags

    return df


# =============================================================================
# Plot
# =============================================================================
def plot_rx_ry(df, centers, out_prefix):

    plt.figure(figsize=(6,6))

    color_map = {
        "PASS": "#4CAF50",
        "OUTLIER": "#FF9800",
        "INCONSISTENT": "#F44336",
        "MOSAIC_CANDIDATE": "#9C27B0"
    }

    for flag in df["cohort_qc_flag"].unique():
        sub = df[df["cohort_qc_flag"] == flag]
        plt.scatter(
            sub["chrX_ratio"],
            sub["chrY_ratio"],
            label=flag,
            alpha=0.7,
            color=color_map.get(flag, "grey")
        )

    # cluster centers
    for cid, (x, y) in centers.items():
        plt.scatter(x, y, marker="X", s=120, color="black")
        plt.text(x, y, f"C{cid}")

    # expected
    for k, (x, y) in EXPECTED.items():
        plt.scatter(x, y, marker="x", s=100)
        plt.text(x, y, k)

    plt.axvline(0.5, linestyle="dotted", color="grey")
    plt.axvline(1.0, linestyle="dotted", color="grey")
    plt.axhline(0.5, linestyle="dotted", color="grey")

    plt.xlabel("RX")
    plt.ylabel("RY")
    plt.title("Cohort clustering (DBSCAN)")

    plt.legend()
    plt.xlim(0, 2)
    plt.ylim(0, 2)

    plt.savefig(f"{out_prefix}_rx_ry.png", dpi=150)
    plt.close()


# =============================================================================
# Main
# =============================================================================
def main():

    args = parse_args()

    log("Loading data")
    df = load_data(args.input_tsvs)

    log("Running clustering")
    df = run_clustering(df)

    log("Computing cluster centers")
    centers = compute_cluster_centers(df)

    log("Assigning QC flags")
    df = assign_qc(df, centers)

    log("Plotting")
    plot_rx_ry(df, centers, args.out_prefix)

    output = df[[
        "individual",
        "chrX_ratio",
        "chrY_ratio",
        "karyotype",
        "cluster_id",
        "cohort_qc_flag",
        "cohort_confidence",
        "mosaic_flag",
        "arm_flags",
        "arm_confidence"
    ]]

    output.to_csv(f"{args.out_prefix}_qc.tsv", sep="\t", index=False)

    log("Done")


if __name__ == "__main__":
    main()