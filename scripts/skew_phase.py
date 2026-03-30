#!/usr/bin/env python3
"""
skew_phase.py — Skew-aware phasing of chrX variants.

Uses X-inactivation skew signal to orient phase blocks consistently
across the entire X chromosome, then outputs a bgzipped VCF.

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


def load_skew(skew_file):
    """
    Load skew file and determine which haplotype is the active X per phase block.
    Uses H1_Xa and H2_Xa read counts to determine active haplotype.
    Only returns phase blocks where there is a clear skew signal.
    Returns dict: PS -> "H1" or "H2"
    """
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

            ps   = parts[ps_idx]
            h1_xa = float(parts[h1_xa_idx])
            h2_xa = float(parts[h2_xa_idx])

            total = h1_xa + h2_xa
            if total == 0:
                continue

            # Only use phase blocks with meaningful read support
            if total < 5:
                continue

            if h1_xa > h2_xa:
                ps_active[ps] = "H1"
            else:
                ps_active[ps] = "H2"

    print(f"[skew_phase] Loaded {len(ps_active)} phase blocks with skew signal", file=sys.stderr)
    return ps_active


def process_vcf(vcf_in, vcf_out_path, ps_active):
    """
    Process VCF:
    - Only modify chrX variants
    - For heterozygous phased variants in a known phase block,
      orient so that the active X haplotype carries the ref (0) allele
      on H1 consistently
    - Write bgzipped output
    """
    flipped   = 0
    kept      = 0
    skipped   = 0
    chrx_total = 0

    with gzip.open(vcf_in, "rt") as fin, gzip.open(vcf_out_path, "wt") as fout:

        for line in fin:

            # Pass through header lines unchanged
            if line.startswith("#"):
                fout.write(line)
                continue

            parts = line.strip().split("\t")
            chrom = parts[0]

            # Only modify chrX variants — pass everything else through unchanged
            if chrom not in ("chrX", "X"):
                fout.write(line)
                continue

            chrx_total += 1

            fmt    = parts[8].split(":")
            sample = parts[9].split(":")

            # Skip if no PS or GT in format
            if "PS" not in fmt or "GT" not in fmt:
                fout.write(line)
                skipped += 1
                continue

            ps_val = sample[fmt.index("PS")]
            gt     = sample[fmt.index("GT")]

            # Skip missing genotypes or phase blocks not in skew data
            if gt in ("./.", ".|.", ".") or ps_val not in ps_active:
                fout.write(line)
                skipped += 1
                continue

            active = ps_active[ps_val]

            # Only flip heterozygous phased genotypes
            # If H2 is active, flip so active haplotype is consistently H1
            if gt in ("0|1", "1|0"):
                if active == "H2":
                    # Flip: swap alleles so H1 = active X
                    new_gt = "1|0" if gt == "0|1" else "0|1"
                    sample[fmt.index("GT")] = new_gt
                    parts[9] = ":".join(sample)
                    flipped += 1
                else:
                    kept += 1
            else:
                # Homozygous or unphased — leave unchanged
                kept += 1

            fout.write("\t".join(parts) + "\n")

    print(f"[skew_phase] chrX variants processed : {chrx_total}", file=sys.stderr)
    print(f"[skew_phase] Variants flipped         : {flipped}",    file=sys.stderr)
    print(f"[skew_phase] Variants kept as-is      : {kept}",       file=sys.stderr)
    print(f"[skew_phase] Variants skipped (no PS) : {skipped}",    file=sys.stderr)


def index_vcf(vcf_path):
    """bgzip is already handled by gzip.open wt — run tabix to index."""
    print(f"[skew_phase] Indexing {vcf_path} with tabix...", file=sys.stderr)
    result = subprocess.run(["tabix", "-p", "vcf", vcf_path], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[skew_phase] WARNING: tabix failed: {result.stderr}", file=sys.stderr)
        print(f"[skew_phase] You can index manually with: tabix -p vcf {vcf_path}", file=sys.stderr)
    else:
        print(f"[skew_phase] Index created: {vcf_path}.tbi", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Skew-aware chrX phasing — orients phase blocks using X-inactivation skew signal"
    )
    parser.add_argument("--vcf",  required=True, help="Input phased VCF (bgzipped .vcf.gz)")
    parser.add_argument("--skew", required=True, help="Skew TSV (gzipped, from SkewX pipeline)")
    parser.add_argument("--out",  required=True, help="Output VCF path (will be bgzipped .vcf.gz)")

    args = parser.parse_args()

    # Ensure output ends in .vcf.gz
    out_path = args.out
    if not out_path.endswith(".vcf.gz"):
        out_path = out_path + ".vcf.gz"
        print(f"[skew_phase] Output path adjusted to: {out_path}", file=sys.stderr)

    print(f"[skew_phase] Loading skew data from: {args.skew}", file=sys.stderr)
    ps_active = load_skew(args.skew)

    print(f"[skew_phase] Processing VCF: {args.vcf}", file=sys.stderr)
    process_vcf(args.vcf, out_path, ps_active)

    # Attempt tabix indexing
    index_vcf(out_path)

    print(f"[skew_phase] Done. Output: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()