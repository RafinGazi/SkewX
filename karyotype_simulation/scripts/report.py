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
EXPECTED_SKIP = {"XY", "XO", "XXX", "XXY", "XYY"}


def generate_report(results_dir):

    tsv_files = sorted(glob.glob(f"{results_dir}/*_karyotype.tsv"))

    if not tsv_files:
        log("No karyotype TSVs found — run infer_karyotype.py first")
        return

    df = pd.concat([pd.read_csv(f, sep="\t") for f in tsv_files], ignore_index=True)

    df["arm_flags"] = df["arm_flags"].fillna("none")
    df["xp_xq_ratio"] = df.get("xp_xq_ratio", pd.Series([float("nan")] * len(df)))

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
            f"<td>{r['status']}</td>"
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
        elif sample in EXPECTED_SKIP:
            match = "✅" if inferred == expected and qc.startswith("skipped") else "❌"
        else:
            match = "✅" if inferred == expected else "❌"

        validation_rows += (
            f"<tr>"
            f"<td>{sample}</td>"
            f"<td>{expected}</td>"
            f"<td>{inferred}</td>"
            f"<td>{arm_flags}</td>"
            f"<td>{qc}</td>"
            f"<td style='text-align:center'>{match}</td>"
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

    <h2>Overview</h2>
    <p>
    This report validates karyotype inference using simulated BAM files with known chromosomal configurations.
    The pipeline estimates chromosome copy number from coverage (RX/RY) and detects structural abnormalities
    using Xp vs Xq arm imbalance.
    </p>

    <h2>Simulation Design</h2>
    <table border="1" cellpadding="6">
    <tr><th>Sample</th><th>Description</th><th>Expected Karyotype</th><th>Expected QC</th></tr>

    <tr><td>XX</td><td>Normal diploid female</td><td>XX</td><td>pass</td></tr>

    <tr><td>XY</td><td>Male (1 X, 1 Y)</td><td>XY</td><td>skipped</td></tr>
    <tr><td>XO</td><td>Monosomy X (single X chromosome)</td><td>XO</td><td>skipped</td></tr>
    <tr><td>XXX</td><td>Trisomy X</td><td>XXX</td><td>skipped</td></tr>
    <tr><td>XXY</td><td>Klinefelter syndrome (2 X, 1 Y)</td><td>XXY</td><td>skipped</td></tr>
    <tr><td>XYY</td><td>Extra Y chromosome</td><td>XYY</td><td>skipped</td></tr>

    <tr><td>Xp_deleted</td><td>Xp arm removed (coverage drops to ~0)</td><td>unknown</td><td>flagged</td></tr>
    <tr><td>Xq_deleted</td><td>Xq arm removed (coverage drops to ~0)</td><td>unknown</td><td>flagged</td></tr>
    <tr><td>X_partial</td><td>Partial Xp reduction (~50% coverage)</td><td>unknown</td><td>flagged</td></tr>

    </table>

    <h2>Expected Behaviour</h2>
    <ul>
    <li><b>XX:</b> Balanced X coverage → passes QC</li>
    <li><b>Non-XX (XY, XO, XXX, XXY, XYY):</b> Identified by RX/RY → skipped</li>
    <li><b>Structural abnormalities:</b> Detected via Xp/Xq imbalance → flagged</li>
    <li><b>Key principle:</b> XO shows uniform reduction across X, whereas deletions show arm-specific drops</li>
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