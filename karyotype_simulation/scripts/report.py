#!/usr/bin/env python3

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


EXPECTED_MAP = {
    "XX": "XX",
    "XY": "XY",
    "XO": "XO",
    "XXX": "XXX",
    "XXY": "XXY",
    "XYY": "XYY",
    "Xp_deleted": "unknown",
    "Xq_deleted": "unknown",
    "X_partial": "unknown",
}

STRUCTURAL = {"Xp_deleted", "Xq_deleted", "X_partial"}
EXPECTED_SKIP = {"XY", "XO", "XYY"}
EXPECTED_PASS = {"XX", "XXX", "XXY"}

def generate_report(results_dir):

    tsv_files = sorted(glob.glob(f"{results_dir}/*_karyotype.tsv"))

    if not tsv_files:
        log("No karyotype TSVs found — run infer_karyotype.py first")
        return

    df = pd.concat([pd.read_csv(f, sep="\t") for f in tsv_files], ignore_index=True)

    df["arm_flags"] = df["arm_flags"].fillna("none")
    df["xp_xq_ratio"] = df.get("xp_xq_ratio", pd.Series([float("nan")] * len(df)))

    # =========================
    # Summary
    # =========================
    total = len(df)
    correct = 0

    for _, r in df.iterrows():
        sample = r["individual"].replace("_sorted", "")
        qc = r["qc_flag"]

        if sample in STRUCTURAL and qc.startswith("flagged"):
            correct += 1
        elif sample in EXPECTED_SKIP and qc.startswith("skipped"):
            correct += 1
        elif sample in ["XX", "XXX", "XXY"] and qc == "pass":
            correct += 1

    summary_html = f"""
    <p style="font-size:16px;">
    <b>Summary:</b> {correct}/{total} samples matched expected validation outcome
    </p>
    """

    # =========================
    # Results table
    # =========================
    results_rows = ""
    for _, r in df.iterrows():
        qc = r["qc_flag"]

        if qc == "pass":
            bg = 'style="background:#e8f5e9;"'
        elif qc.startswith("skipped"):
            bg = 'style="background:#fff3e0;"'
        elif qc.startswith("flagged"):
            bg = 'style="background:#fce4ec;"'
        else:
            bg = ""

        results_rows += (
            f"<tr {bg}>"
            f"<td>{r['individual']}</td>"
            f"<td>{r['karyotype']}</td>"
            f"<td>{r.get('status', 'n/a')}</td>"
            f"<td>{r['chrX_ratio']:.3f}</td>"
            f"<td>{r['chrY_ratio']:.3f}</td>"
            f"<td>{r['auto_mean_cov']:.1f}x</td>"
            f"<td>{r['arm_flags']}</td>"
            f"<td>{qc}</td>"
            f"</tr>"
        )

    # =========================
    # Validation table
    # =========================
    validation_rows = ""
    for _, r in df.iterrows():
        sample = r["individual"].replace("_sorted", "")
        expected = EXPECTED_MAP.get(sample, "unknown")
        inferred = r["karyotype"]
        arm_flags = r["arm_flags"]
        qc = r["qc_flag"]

        if sample in STRUCTURAL:
            match = "✅" if qc.startswith("flagged") else "❌"
        elif sample in EXPECTED_PASS:
            match = "✅" if qc == "pass" else "❌"
        elif sample in EXPECTED_SKIP:
            match = "✅" if inferred == expected and qc.startswith("skipped") else "❌"
        else:
            match = "❌"

        color = "#e8f5e9" if match == "✅" else "#ffebee"
        validation_rows += (
            f"<tr>"
            f"<td>{sample}</td>"
            f"<td>{expected}</td>"
            f"<td>{inferred}</td>"
            f"<td>{arm_flags}</td>"
            f"<td>{qc}</td>"
            f"<td style='text-align:center; background:{color}'>{match}</td>"
            f"</tr>"
        )

    # =========================
    # Per-sample plots
    # =========================
    per_sample_plots = ""

    for _, r in df.iterrows():
        name = r["individual"]
        karyotype = r["karyotype"]
        qc = r["qc_flag"]
        arm_flags = r["arm_flags"]

        xp_xq_str = (
            f"{r['xp_xq_ratio']:.3f}"
            if pd.notna(r.get("xp_xq_ratio"))
            else "n/a"
        )

        if qc.startswith("skipped"):
            notice = f"<p><b>⚠️ Skipped:</b> {karyotype}</p>"
        elif qc.startswith("flagged"):
            notice = f"<p><b>🚩 Structural:</b> {arm_flags}</p>"
        else:
            notice = "<p><b>✅ Passed</b></p>"

        coverage_plot = img_tag(f"{results_dir}/{name}_coverage_ratios.png")
        rxry_plot = img_tag(f"{results_dir}/{name}_rx_ry.png")
        genome_plot = img_tag(f"{results_dir}/{name}_genome_coverage.png")
        xp_xq_plot = img_tag(f"{results_dir}/{name}_xp_xq.png")

        per_sample_plots += (
            f"<h3>{name} — {karyotype}</h3>"
            f"{notice}"
            f"<p>RX={r['chrX_ratio']:.3f} | RY={r['chrY_ratio']:.3f} | "
            f"Auto={r['auto_mean_cov']:.1f} | Xp/Xq={xp_xq_str}</p>"
            f"{coverage_plot}{rxry_plot}{xp_xq_plot}{genome_plot}"
        )

    html = f"""
    <html>
    <body style="font-family: Arial; margin: 40px; max-width: 1400px;">

    <h1>Karyotype Validation Report</h1>
    {summary_html}

    <h2>1. Source Data</h2>
    <table border="1" cellpadding="6">
    <tr><th>Component</th><th>Source</th></tr>
    <tr><td>Autosomes</td><td>GM19312 chr18 + chr21</td></tr>
    <tr><td>chrX</td><td>GM19312 chrX</td></tr>
    <tr><td>chrY</td><td>GM19312 chrY</td></tr>
    <tr><td>Structural simulations</td><td>GM19462 chrX (Xp/Xq deletions)</td></tr>
    </table>

    <br>

    <h2>2. Simulation Design</h2>
    <table border="1" cellpadding="6">
    <tr>
    <th>Sample</th>
    <th>Composition</th>
    <th>Expected inference</th>
    <th>Expected pipeline action</th>
    </tr>

    <tr><td>XX</td><td>AUTO + X + X</td><td>XX</td><td>Pass</td></tr>
    <tr><td>XY</td><td>AUTO + X + Y</td><td>XY</td><td>Skip</td></tr>
    <tr><td>XO</td><td>AUTO + 0.5X</td><td>XO</td><td>Skip</td></tr>
    <tr><td>XXX</td><td>AUTO + X + X + X</td><td>XXX</td><td>Skip</td></tr>
    <tr><td>XXY</td><td>AUTO + X + X + Y</td><td>XXY</td><td>Skip</td></tr>
    <tr><td>XYY</td><td>AUTO + X + Y + Y</td><td>XYY</td><td>Skip</td></tr>

    <tr><td>Xp_deleted</td><td>AUTO + Xq only (Xp removed)</td><td>unknown + Xp_deletion</td><td>Flag</td></tr>
    <tr><td>Xq_deleted</td><td>AUTO + Xp only (Xq removed)</td><td>unknown + Xq_deletion</td><td>Flag</td></tr>
    <tr><td>X_partial</td><td>AUTO + Xp@50% + Xq full</td><td>unknown + Xp_partial_deletion</td><td>Flag</td></tr>
    </table>

    <br>

    <h2>3. Interpretation Logic</h2>
    <ul>
    <li><b>Karyotype inference:</b> Based on relative coverage (RX and RY)</li>
    <li><b>XX:</b> Balanced X coverage → passes QC</li>
    <li><b>XXX / XXY:</b> X-inactivation present → included in analysis</li>
    <li><b>XY / XO / XYY:</b> insufficient X dosage → skipped</li>
    <li><b>Structural abnormalities:</b> detected using Xp vs Xq imbalance</li>
    <li><b>Structural output:</b> Reported as "unknown" karyotype with arm-specific flags</li>
    <li><b>Key distinction:</b>
        <ul>
            <li>XO → uniformly reduced X coverage</li>
            <li>Deletions → arm-specific drop (Xp or Xq)</li>
        </ul>
    </li>
    </ul>

    <h2>Results</h2>
    <table border="1">
        <tr>
            <th>Sample</th><th>Karyotype</th><th>Status</th>
            <th>RX</th><th>RY</th><th>Auto</th><th>Flags</th><th>QC</th>
        </tr>
        {results_rows}
    </table>

    <h2>Validation</h2>
    <table border="1">
        <tr>
            <th>Sample</th><th>Expected</th><th>Observed</th>
            <th>Flags</th><th>QC</th><th>Match</th>
        </tr>
        {validation_rows}
    </table>

    <h2>Per Sample</h2>
    {per_sample_plots}

    </body>
    </html>
    """

    out_path = f"{results_dir}/report.html"
    with open(out_path, "w") as f:
        f.write(html)

    log(f"Created {out_path}")


def main():
    generate_report("../results")


if __name__ == "__main__":
    main()