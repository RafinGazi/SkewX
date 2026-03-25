#!/usr/bin/env python3

import pandas as pd
import glob
import os

def log(msg):
    print(f"[report] {msg}")

# =========================
# Individual reports (unchanged)
# =========================
def generate_individual_reports(results_dir):

    files = glob.glob(f"{results_dir}/*_karyotype.tsv")

    for f in files:
        df = pd.read_csv(f, sep="\t")
        row = df.iloc[0]

        sample = row["individual"]

        html = f"""
        <html>
        <body style="font-family: Arial; margin: 40px;">

        <h1>Karyotype Report: {sample}</h1>

        <h2>Summary</h2>
        <p><b>Karyotype:</b> {row['karyotype']}</p>
        <p><b>Confidence:</b> {row['karyotype_conf']}</p>

        <h2>Copy Numbers</h2>
        <p>X: {row['x_copies']} | Y: {row['y_copies']}</p>

        <h2>Coverage Ratios</h2>
        <p>RX: {row['chrX_ratio']}</p>
        <p>RY: {row['chrY_ratio']}</p>

        <h2>Structural Signals</h2>
        <p>Arm flags: {row['arm_flags']}</p>
        <p>Mosaic: {row['mosaic_status']}</p>

        <h2>Plots</h2>
        <img src="{sample}_coverage_ratios.png" width="400"><br><br>
        <img src="{sample}_chrX_coverage.png" width="600"><br><br>
        <img src="{sample}_rx_ry.png" width="400">

        </body>
        </html>
        """

        out_path = f"{results_dir}/{sample}_report.html"

        with open(out_path, "w") as out:
            out.write(html)

        log(f"Created {out_path}")


# =========================
# Cohort report (ENHANCED)
# =========================
def generate_cohort_report(results_dir):

    cohort_file = f"{results_dir}/cohort_qc.tsv"

    if not os.path.exists(cohort_file):
        log("No cohort_qc.tsv found — skipping cohort report")
        return

    df = pd.read_csv(cohort_file, sep="\t")

    # =========================
    # Expected mapping
    # =========================
    expected_map = {
        "XX": "XX",
        "XY": "XY",
        "XO": "XO",
        "XXX": "XXX",
        "XXY": "XXY",
        "XYY": "XYY"
    }

    validation_rows = ""
    for _, r in df.iterrows():

        sample = r["individual"].replace("_sorted", "")
        expected = expected_map.get(sample, "unknown")
        inferred = r["karyotype"]

        match = "✅" if expected == inferred else "❌"

        validation_rows += f"""
        <tr>
            <td>{sample}</td>
            <td>{expected}</td>
            <td>{inferred}</td>
            <td>{r['cohort_qc_flag']}</td>
            <td>{match}</td>
        </tr>
        """

    # =========================
    # Main QC table
    # =========================
    rows = ""
    for _, r in df.iterrows():
        rows += f"""
        <tr>
            <td>{r['individual']}</td>
            <td>{r['karyotype']}</td>
            <td>{r['chrX_ratio']:.3f}</td>
            <td>{r['chrY_ratio']:.3f}</td>
            <td>{r['cohort_qc_flag']}</td>
            <td>{r['cohort_confidence']}</td>
        </tr>
        """

    html = f"""
    <html>
    <body style="font-family: Arial; margin: 40px;">

    <h1>Karyotype Simulation & Validation Report</h1>

    <h2>1. Source Data</h2>
    <table border="1" cellpadding="6">
        <tr><th>Component</th><th>Source</th></tr>
        <tr><td>Autosomes</td><td>GM19312 chr18 + chr21</td></tr>
        <tr><td>chrX</td><td>GM19312 chrX</td></tr>
        <tr><td>chrY</td><td>GM19312 chrY</td></tr>
    </table>

    <h2>2. Simulation Design</h2>
    <table border="1" cellpadding="6">
        <tr><th>Simulated Sample</th><th>Composition</th></tr>
        <tr><td>XX</td><td>AUTO + X + X</td></tr>
        <tr><td>XY</td><td>AUTO + X + Y</td></tr>
        <tr><td>XO</td><td>AUTO + 0.5 X</td></tr>
        <tr><td>XXX</td><td>AUTO + X + X + X</td></tr>
        <tr><td>XXY</td><td>AUTO + X + X + Y</td></tr>
        <tr><td>XYY</td><td>AUTO + X + Y + Y</td></tr>
    </table>

    <h2>3. Cohort QC Results</h2>
    <table border="1" cellpadding="6">
        <tr>
            <th>Sample</th>
            <th>Karyotype</th>
            <th>RX</th>
            <th>RY</th>
            <th>QC Flag</th>
            <th>Confidence</th>
        </tr>
        {rows}
    </table>

    <h2>4. Validation (Expected vs Observed)</h2>
    <table border="1" cellpadding="6">
        <tr>
            <th>Sample</th>
            <th>Expected</th>
            <th>Inferred</th>
            <th>QC</th>
            <th>Match</th>
        </tr>
        {validation_rows}
    </table>

    <h2>Clustering Plot</h2>
    <img src="cohort_rx_ry.png" width="600">

    </body>
    </html>
    """

    out_path = f"{results_dir}/cohort_report.html"

    with open(out_path, "w") as out:
        out.write(html)

    log(f"Created {out_path}")


# =========================
# Main
# =========================
def main():

    results_dir = "../results"

    log("Generating individual reports...")
    generate_individual_reports(results_dir)

    log("Generating cohort report...")
    generate_cohort_report(results_dir)

    log("All reports generated successfully.")


if __name__ == "__main__":
    main()