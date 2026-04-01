#!/usr/bin/env python3

# =============================================================================
# report.py
#
# Validation report for karyotype simulation.
# Reads per-sample TSVs directly — no cohort file needed.
#
# For each sample shows:
#   - Karyotype call + QC flag
#   - Whether it passed or was skipped and why
#   - All plots as proof
# =============================================================================

import os
import glob
import base64
import pandas as pd

def log(msg):
    print(f"[report] {msg}")


def img_tag(path, width=600):
    if not os.path.exists(path):
        return f'<p style="color:grey;">[plot not found: {os.path.basename(path)}]</p>'
    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode("utf-8")
    return f'<img src="data:image/png;base64,{data}" width="{width}" style="display:block;margin:8px 0;">'


# Expected karyotype for each simulated sample
EXPECTED_MAP = {
    "XX":         "XX",
    "XY":         "XY",
    "XO":         "XO",
    "XXX":        "XXX",
    "XXY":        "XXY",
    "XYY":        "XYY",
    "Xp_deleted": "XX",   # structurally abnormal XX
    "Xq_deleted": "XX",   # structurally abnormal XX
    "X_partial":  "XX",   # structurally abnormal XX
}

# Structural samples — karyotype should be XX AND arm_flags should fire
STRUCTURAL = {"Xp_deleted", "Xq_deleted", "X_partial"}

# Samples that should be skipped by the pipeline (not XX)
EXPECTED_SKIP = {"XY", "XO", "XXX", "XXY", "XYY"}


def generate_report(results_dir):

    tsv_files = sorted(glob.glob(f"{results_dir}/*_karyotype.tsv"))

    if not tsv_files:
        log("No karyotype TSVs found — run infer_karyotype.py first")
        return

    df = pd.concat([pd.read_csv(f, sep="\t") for f in tsv_files], ignore_index=True)
    df["arm_flags"]   = df["arm_flags"].fillna("none")
    df["xp_xq_ratio"] = df.get("xp_xq_ratio", pd.Series([float("nan")] * len(df)))

    # =========================
    # Results table
    # =========================
    results_rows = ""
    for _, r in df.iterrows():
        qc      = r["qc_flag"]
        bg      = ""
        if qc == "pass":
            bg = 'style="background:#e8f5e9;"'
        elif qc.startswith("skipped"):
            bg = 'style="background:#fff3e0;"'
        elif qc.startswith("flagged"):
            bg = 'style="background:#fce4ec;"'

        results_rows += f"""
        <tr {bg}>
            <td>{r['individual']}</td>
            <td>{r['karyotype']}</td>
            <td>{r['chrX_ratio']:.3f}</td>
            <td>{r['chrY_ratio']:.3f}</td>
            <td>{r['auto_mean_cov']:.1f}x</td>
            <td>{r['arm_flags']}</td>
            <td>{qc}</td>
        </tr>
        """

    # =========================
    # Validation table
    # =========================
    validation_rows = ""
    for _, r in df.iterrows():
        sample    = r["individual"].replace("_sorted", "")
        expected  = EXPECTED_MAP.get(sample, "unknown")
        inferred  = r["karyotype"]
        arm_flags = r["arm_flags"]
        qc        = r["qc_flag"]

        if sample in STRUCTURAL:
            # Must be called XX and arm_flags must have fired
            match = "✅" if inferred == "XX" and arm_flags != "none" else "❌"
        elif sample in EXPECTED_SKIP:
            # Should be skipped — check karyotype matches and qc_flag starts with skipped
            match = "✅" if inferred == expected and qc.startswith("skipped") else "❌"
        else:
            match = "✅" if inferred == expected else "❌"

        validation_rows += f"""
        <tr>
            <td>{sample}</td>
            <td>{expected}</td>
            <td>{inferred}</td>
            <td>{arm_flags}</td>
            <td>{qc}</td>
            <td style="font-size:16px; text-align:center;">{match}</td>
        </tr>
        """

    # =========================
    # Per-sample plots
    # =========================
    per_sample_plots = ""
    for _, r in df.iterrows():
        name      = r["individual"]
        sample    = name.replace("_sorted", "")
        karyotype = r["karyotype"]
        qc        = r["qc_flag"]
        arm_flags = r["arm_flags"]
        xp_xq_str = f"{r['xp_xq_ratio']:.3f}" if pd.notna(r.get("xp_xq_ratio")) else "n/a"

        # Skip notice for non-XX samples
        if qc.startswith("skipped"):
            skip_notice = f"""
            <div style="background:#fff3e0; border-left:4px solid #ff9800;
                        padding:10px 16px; margin-bottom:12px; border-radius:4px;">
                <b>⚠️ Sample skipped by pipeline</b><br>
                Inferred karyotype: <b>{karyotype}</b> — not XX.<br>
                Skipped at karyotype check stage. Plots below are provided as proof.
            </div>
            """
        elif qc.startswith("flagged"):
            skip_notice = f"""
            <div style="background:#fce4ec; border-left:4px solid #e91e63;
                        padding:10px 16px; margin-bottom:12px; border-radius:4px;">
                <b>🚩 Structural flag detected</b><br>
                Karyotype: <b>{karyotype}</b> — arm flags: <b>{arm_flags}</b>.<br>
                Sample flagged for structural abnormality.
            </div>
            """
        else:
            skip_notice = f"""
            <div style="background:#e8f5e9; border-left:4px solid #4caf50;
                        padding:10px 16px; margin-bottom:12px; border-radius:4px;">
                <b>✅ Sample passed karyotype check</b><br>
                Karyotype: <b>{karyotype}</b> — proceeding to X-inactivation analysis.
            </div>
            """

        coverage_plot = img_tag(f"{results_dir}/{name}_coverage_ratios.png", width=500)
        rxry_plot     = img_tag(f"{results_dir}/{name}_rx_ry.png", width=400)
        genome_plot   = img_tag(f"{results_dir}/{name}_genome_coverage.png", width=900)
        xp_xq_plot    = img_tag(f"{results_dir}/{name}_xp_xq.png", width=400)

        per_sample_plots += f"""
        <div style="border:1px solid #ccc; padding:16px; margin-bottom:24px; border-radius:6px;">
            <h3>{name} — {karyotype}</h3>
            {skip_notice}
            <p>
                <b>RX:</b> {r['chrX_ratio']:.3f} &nbsp;|&nbsp;
                <b>RY:</b> {r['chrY_ratio']:.3f} &nbsp;|&nbsp;
                <b>Autosome mean:</b> {r['auto_mean_cov']:.1f}x &nbsp;|&nbsp;
                <b>Xp/Xq ratio:</b> {xp_xq_str} &nbsp;|&nbsp;
                <b>Arm flags:</b> {arm_flags} &nbsp;|&nbsp;
                <b>QC:</b> {qc}
            </p>
            <div style="display:flex; gap:24px; flex-wrap:wrap;">
                <div>{coverage_plot}<p style="font-size:12px;color:grey;">Coverage ratios</p></div>
                <div>{rxry_plot}<p style="font-size:12px;color:grey;">RX vs RY</p></div>
                <div>{xp_xq_plot}<p style="font-size:12px;color:grey;">Xp vs Xq</p></div>
            </div>
            {genome_plot}
            <p style="font-size:12px;color:grey;">Genome coverage profile (chr18, chr21, chrX, chrY)</p>
        </div>
        """

    # =========================
    # HTML
    # =========================
    html = f"""
    <html>
    <body style="font-family: Arial; margin: 40px; max-width: 1400px;">

    <h1>Karyotype Simulation Validation Report</h1>

    <h2>1. Source Data</h2>
    <table border="1" cellpadding="6">
        <tr><th>Component</th><th>Source</th></tr>
        <tr><td>Autosomes</td><td>GM19312 chr18 + chr21</td></tr>
        <tr><td>chrX</td><td>GM19312 chrX</td></tr>
        <tr><td>chrY</td><td>GM19312 chrY</td></tr>
        <tr><td>Structural sims</td><td>GM19462 chrX (Xp/Xq deletions)</td></tr>
    </table>

    <h2>2. Simulation Design</h2>
    <table border="1" cellpadding="6">
        <tr><th>Sample</th><th>Composition</th><th>Expected karyotype call</th><th>Expected pipeline action</th></tr>
        <tr><td>XX</td><td>AUTO + X + X</td><td>XX</td><td>Pass</td></tr>
        <tr><td>XY</td><td>AUTO + X + Y</td><td>XY</td><td>Skip</td></tr>
        <tr><td>XO</td><td>AUTO + 0.5X</td><td>XO</td><td>Skip</td></tr>
        <tr><td>XXX</td><td>AUTO + X + X + X</td><td>XXX</td><td>Skip</td></tr>
        <tr><td>XXY</td><td>AUTO + X + X + Y</td><td>XXY</td><td>Skip</td></tr>
        <tr><td>XYY</td><td>AUTO + X + Y + Y</td><td>XYY</td><td>Skip</td></tr>
        <tr><td>Xp_deleted</td><td>AUTO + Xq only (Xp removed)</td><td>XX + Xp_deletion flag</td><td>Flag</td></tr>
        <tr><td>Xq_deleted</td><td>AUTO + Xp only (Xq removed)</td><td>XX + Xq_deletion flag</td><td>Flag</td></tr>
        <tr><td>X_partial</td><td>AUTO + Xp@50% + Xq full</td><td>XX + Xp_partial_deletion flag</td><td>Flag</td></tr>
    </table>

    <h2>3. Results</h2>
    <table border="1" cellpadding="6">
        <tr>
            <th>Sample</th><th>Karyotype</th><th>RX</th><th>RY</th>
            <th>Autosome Mean</th><th>Arm Flags</th><th>QC Flag</th>
        </tr>
        {results_rows}
    </table>

    <h2>4. Validation (Expected vs Observed)</h2>
    <table border="1" cellpadding="6">
        <tr>
            <th>Sample</th><th>Expected</th><th>Inferred</th>
            <th>Arm Flags</th><th>QC</th><th>Match</th>
        </tr>
        {validation_rows}
    </table>

    <h2>5. Per-Sample Plots</h2>
    {per_sample_plots}

    </body>
    </html>
    """

    out_path = f"{results_dir}/report.html"
    with open(out_path, "w") as f:
        f.write(html)

    log(f"Created {out_path}")


def main():
    results_dir = "../results"
    generate_report(results_dir)


if __name__ == "__main__":
    main()