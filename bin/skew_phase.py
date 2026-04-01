#!/usr/bin/env python3
"""
skew_phase.py — Skew-aware phasing of chrX variants.

"""

import gzip
import argparse
import sys


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

    print(f"[skew_phase] Loaded {len(ps_active)} phase blocks", file=sys.stderr)
    return ps_active


def process_vcf(vcf_in, vcf_out, ps_active):
    flipped = kept = skipped = chrx_total = 0

    with gzip.open(vcf_in, "rt") as fin, open(vcf_out, "w") as fout:
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

    print(f"[skew_phase] chrX processed : {chrx_total}", file=sys.stderr)
    print(f"[skew_phase] flipped        : {flipped}", file=sys.stderr)
    print(f"[skew_phase] kept           : {kept}", file=sys.stderr)
    print(f"[skew_phase] skipped        : {skipped}", file=sys.stderr)

    total = flipped + kept

    flip_rate = (flipped / total) if total > 0 else 0

    with open(vcf_out + ".metrics.txt", "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"flipped\t{flipped}\n")
        f.write(f"kept\t{kept}\n")
        f.write(f"skipped\t{skipped}\n")
        f.write(f"flip_rate\t{flip_rate:.4f}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Skew-aware chrX phasing"
    )

    parser.add_argument("--vcf", required=True)
    parser.add_argument("--skew", required=True)
    parser.add_argument("--out", required=True)

    args = parser.parse_args()

    print(f"[skew_phase] Loading skew...", file=sys.stderr)
    ps_active = load_skew(args.skew)

    print(f"[skew_phase] Processing VCF...", file=sys.stderr)
    process_vcf(args.vcf, args.out, ps_active)

    print(f"[skew_phase] Done → {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()