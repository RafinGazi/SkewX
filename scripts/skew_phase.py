#!/usr/bin/env python3
"""
skew_phase.py — Skew-aware phasing of chrX variants.

Uses X-inactivation skew signal to orient phase blocks consistently
across the entire X chromosome, then outputs a bgzipped + tabix indexed VCF.

Usage:
    python3 skew_phase.py --vcf input.phased.vcf.gz \
                          --skew sample_skew.tsv.gz \
                          --out output.skew_phased.vcf.gz
"""

import gzip
import argparse
import subprocess
import sys
import os
import shutil


def find_tool(name):
    path = shutil.which(name)
    if path:
        return path
    fallbacks = [
        "/opt/apps/eb/software/CellRanger/3.1.0/miniconda-cr-cs/4.3.21-miniconda-cr-cs-c10/bin",
    ]
    for d in fallbacks:
        candidate = os.path.join(d, name)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return None


def load_skew(skew_file):
    ps_active = {}
    with gzip.open(skew_file, "rt") as f:
        header = next(f).strip().split("\t")
        h1_xa_idx = header.index("H1_Xa")
        h2_xa_idx = header.index("H2_Xa")
        ps_idx    = header.index("PS")
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) <= max(h1_xa_idx, h2_xa_idx, ps_idx):
                continue
            ps    = parts[ps_idx]
            h1_xa = float(parts[h1_xa_idx])
            h2_xa = float(parts[h2_xa_idx])
            total = h1_xa + h2_xa
            if total < 5:
                continue
            ps_active[ps] = "H1" if h1_xa > h2_xa else "H2"
    print(f"[skew_phase] Loaded {len(ps_active)} phase blocks with skew signal", file=sys.stderr)
    return ps_active


def process_vcf(vcf_in, vcf_out_plain, ps_active):
    flipped = kept = skipped = chrx_total = 0
    with gzip.open(vcf_in, "rt") as fin, open(vcf_out_plain, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            parts = line.strip().split("\t")
            chrom = parts[0]
            if chrom not in ("chrX", "X"):
                fout.write(line)
                continue
            chrx_total += 1
            fmt    = parts[8].split(":")
            sample = parts[9].split(":")
            if "PS" not in fmt or "GT" not in fmt:
                fout.write(line)
                skipped += 1
                continue
            ps_val = sample[fmt.index("PS")]
            gt     = sample[fmt.index("GT")]
            if gt in ("./.", ".|.", ".") or ps_val not in ps_active:
                fout.write(line)
                skipped += 1
                continue
            active = ps_active[ps_val]
            if gt in ("0|1", "1|0"):
                if active == "H2":
                    new_gt = "1|0" if gt == "0|1" else "0|1"
                    sample[fmt.index("GT")] = new_gt
                    parts[9] = ":".join(sample)
                    flipped += 1
                else:
                    kept += 1
            else:
                kept += 1
            fout.write("\t".join(parts) + "\n")
    print(f"[skew_phase] chrX variants processed : {chrx_total}", file=sys.stderr)
    print(f"[skew_phase] Variants flipped         : {flipped}",    file=sys.stderr)
    print(f"[skew_phase] Variants kept as-is      : {kept}",       file=sys.stderr)
    print(f"[skew_phase] Variants skipped (no PS) : {skipped}",    file=sys.stderr)


def bgzip_and_index(plain_vcf, bgzip_bin, tabix_bin):
    gz_path = plain_vcf + ".gz"
    print(f"[skew_phase] bgzipping {plain_vcf}...", file=sys.stderr)
    r = subprocess.run([bgzip_bin, "-f", plain_vcf], capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[skew_phase] ERROR: bgzip failed: {r.stderr}", file=sys.stderr)
        sys.exit(1)
    print(f"[skew_phase] Indexing with tabix...", file=sys.stderr)
    r = subprocess.run([tabix_bin, "-p", "vcf", gz_path], capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[skew_phase] ERROR: tabix failed: {r.stderr}", file=sys.stderr)
        sys.exit(1)
    print(f"[skew_phase] Index created: {gz_path}.tbi", file=sys.stderr)
    return gz_path


def main():
    parser = argparse.ArgumentParser(
        description="Skew-aware chrX phasing using X-inactivation skew signal"
    )
    parser.add_argument("--vcf",  required=True, help="Input phased VCF (bgzipped .vcf.gz)")
    parser.add_argument("--skew", required=True, help="Skew TSV (gzipped, from SkewX pipeline)")
    parser.add_argument("--out",  required=True, help="Output VCF path (.vcf.gz)")
    args = parser.parse_args()

    bgzip_bin = find_tool("bgzip")
    tabix_bin = find_tool("tabix")

    if not bgzip_bin:
        print("[skew_phase] ERROR: bgzip not found.", file=sys.stderr)
        sys.exit(1)
    if not tabix_bin:
        print("[skew_phase] ERROR: tabix not found.", file=sys.stderr)
        sys.exit(1)

    print(f"[skew_phase] Using bgzip: {bgzip_bin}", file=sys.stderr)
    print(f"[skew_phase] Using tabix: {tabix_bin}", file=sys.stderr)

    out_path = args.out
    if out_path.endswith(".vcf.gz"):
        plain_vcf = out_path[:-3]
    elif out_path.endswith(".vcf"):
        plain_vcf = out_path
        out_path  = out_path + ".gz"
    else:
        plain_vcf = out_path + ".vcf"
        out_path  = plain_vcf + ".gz"
        print(f"[skew_phase] Output path adjusted to: {out_path}", file=sys.stderr)

    print(f"[skew_phase] Loading skew data from: {args.skew}", file=sys.stderr)
    ps_active = load_skew(args.skew)

    print(f"[skew_phase] Processing VCF: {args.vcf}", file=sys.stderr)
    process_vcf(args.vcf, plain_vcf, ps_active)

    final_path = bgzip_and_index(plain_vcf, bgzip_bin, tabix_bin)

    print(f"[skew_phase] Done. Output : {final_path}", file=sys.stderr)
    print(f"[skew_phase]       Index  : {final_path}.tbi", file=sys.stderr)


if __name__ == "__main__":
    main()