#!/usr/bin/env python3

import pandas as pd
import os
import base64

def log(msg):
    print(f"[report] {msg}")


def img_tag(path, width=600):
    """Embed image as base64 so the HTML is self-contained."""
    if not os.path.exists(path):
        return f'<p style="color:grey;">[plot not found: {os.path.basename(path)}]</p>'
    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode("utf-8")
    return f'<img src="data:image/png;base64,{data}" width="{width}" style="display:block;margin:8px 0;">'


def generate_full_report(results_dir):

    cohort_file = f"{results_dir}/cohort_qc.tsv"

    if not os.path.exists(cohort_file):
        log("No cohort_qc.tsv found — skipping report")
        return

    df = pd.read_csv(cohort_file, sep="\t")

    # =========================
    # Cohort table
    # =========================
    cohort_rows = ""
    for _, r in df.iterrows():
        cohort_rows += f"""
        <tr>
            <td>{r['individual']}</td>
            <td>{r['karyotype']}</td>
            <td>{r['chrX_ratio']:.3f}</td>
            <td>{r['chrY_ratio']:.3f}</td>
            <td>{r.get('arm_flags', 'none')}</td>
            <td>{r['cohort_qc_flag']}</td>
            <td>{r['cohort_confidence']}</td>
        </tr>
        """

    # =========================
    # Validation table
    # =========================
    expected_map = {
        "XX":         "XX",
        "XY":         "XY",
        "XO":         "XO",
        "XXX":        "XXX",
        "XXY":        "XXY",
        "XYY":        "XYY",
        "Xp_deleted": "unknown",
        "Xq_deleted": "unknown",
        "X_partial":  "unknown",
    }

    validation_rows = ""
    for _, r in df.iterrows():
        sample   = r["individual"].replace("_sorted", "")
        expected = expected_map.get(sample, "unknown")
        inferred = r["karyotype"]
        arm_flags = r.get("arm_flags", "none")

        # Match logic: structural samples match if arm_flags fired
        if expected == "unknown":
            match = "✅" if arm_flags not in ("none", "", float("nan")) and str(arm_flags) != "nan" else "❌"
        else:
            match = "✅" if expected == inferred else "❌"

        validation_rows += f"""
        <tr>
            <td>{sample}</td>
            <td>{expected}</td>
            <td>{inferred}</td>
            <td>{arm_flags}</td>
            <td>{r['cohort_qc_flag']}</td>
            <td>{match}</td>
        </tr>
        """

    # =========================
    # Per-sample plots
    # =========================
    per_sample_plots = ""
    for _, r in df.iterrows():
        name = r["individual"]
        coverage_plot = img_tag(f"{results_dir}/{name}_coverage_ratios.png", width=500)
        rxry_plot     = img_tag(f"{results_dir}/{name}_rx_ry.png", width=400)
        chrx_plot     = img_tag(f"{results_dir}/{name}_chrX_coverage.png", width=700)
        xp_xq_plot = img_tag(f"{results_dir}/{name}_xp_xq.png", width=400)

        per_sample_plots += f"""
        <div style="border:1px solid #ccc; padding:16px; margin-bottom:24px; border-radius:6px;">
            <h3>{name} — {r['karyotype']}</h3>
            <p><b>RX:</b> {r['chrX_ratio']:.3f} &nbsp;|&nbsp;
               <b>RY:</b> {r['chrY_ratio']:.3f} &nbsp;|&nbsp;
               <b>Arm flags:</b> {r.get('arm_flags', 'none')} &nbsp;|&nbsp;
               <b>QC:</b> {r['cohort_qc_flag']}</p>
            <div style="display:flex; gap:24px; flex-wrap:wrap;">
                <div>{coverage_plot}<p style="font-size:12px;color:grey;">Coverage ratios</p></div>
                <div>{rxry_plot}<p style="font-size:12px;color:grey;">RX vs RY</p></div>
                <div>{xp_xq_plot}<p style="font-size:12px;color:grey;">Xp vs Xq</p></div>
            </div>
            {chrx_plot}
            <p style="font-size:12px;color:grey;">chrX coverage profile</p>
        </div>
        """

    # =========================
    # HTML
    # =========================
    html = f"""
    <html>
    <body style="font-family: Arial; margin: 40px; max-width: 1200px;">

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
        <tr><th>Sample</th><th>Composition</th></tr>
        <tr><td>XX</td><td>AUTO + X + X</td></tr>
        <tr><td>XY</td><td>AUTO + X + Y</td></tr>
        <tr><td>XO</td><td>AUTO + 0.5X</td></tr>
        <tr><td>XXX</td><td>AUTO + X + X + X</td></tr>
        <tr><td>XXY</td><td>AUTO + X + X + Y</td></tr>
        <tr><td>XYY</td><td>AUTO + X + Y + Y</td></tr>
        <tr><td>Xp_deleted</td><td>AUTO + Xq only (Xp removed)</td></tr>
        <tr><td>Xq_deleted</td><td>AUTO + Xp only (Xq removed)</td></tr>
        <tr><td>X_partial</td><td>AUTO + Xp@50% + Xq full</td></tr>
    </table>

    <h2>3. Cohort QC Results</h2>
    <table border="1" cellpadding="6">
        <tr>
            <th>Sample</th>
            <th>Karyotype</th>
            <th>RX</th>
            <th>RY</th>
            <th>Arm Flags</th>
            <th>QC Flag</th>
            <th>Confidence</th>
        </tr>
        {cohort_rows}
    </table>

    <h2>4. Validation (Expected vs Observed)</h2>
    <table border="1" cellpadding="6">
        <tr>
            <th>Sample</th>
            <th>Expected</th>
            <th>Inferred</th>
            <th>Arm Flags</th>
            <th>QC</th>
            <th>Match</th>
        </tr>
        {validation_rows}
    </table>

    <h2>5. Cohort Clustering</h2>
    {img_tag(f"{results_dir}/cohort_rx_ry.png", width=600)}

    <h2>6. Per-Sample Plots</h2>
    {per_sample_plots}

    </body>
    </html>
    """

    out_path = f"{results_dir}/report.html"
    with open(out_path, "w") as out:
        out.write(html)

    log(f"Created {out_path}")


def main():
    results_dir = "../results"
    generate_full_report(results_dir)


if __name__ == "__main__":
    main()