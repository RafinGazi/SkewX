#!/usr/bin/env python3

import pandas as pd
import glob
import os

def log(msg):
    print(f"[report] {msg}")

# =========================
# Individual reports
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
# Cohort report
# =========================
def generate_cohort_report(results_dir):

    cohort_file = f"{results_dir}/cohort_qc.tsv"

    if not os.path.exists(cohort_file):
        log("No cohort_qc.tsv found — skipping cohort report")
        return

    df = pd.read_csv(cohort_file, sep="\t")

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

    <h1>Cohort QC Report</h1>

    <table border="1" cellpadding="6" cellspacing="0">
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
