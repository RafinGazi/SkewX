import gzip
import argparse

def load_skew(skew_file):
    ps_active = {}

    with gzip.open(skew_file, "rt") as f:
        header = next(f).strip().split("\t")

        h1_xa_idx = header.index("H1_Xa")
        h2_xa_idx = header.index("H2_Xa")
        ps_idx = header.index("PS")

        for line in f:
            parts = line.strip().split("\t")

            ps = parts[ps_idx]
            h1_xa = float(parts[h1_xa_idx])
            h2_xa = float(parts[h2_xa_idx])

            if h1_xa > h2_xa:
                ps_active[ps] = "H1"
            else:
                ps_active[ps] = "H2"

    return ps_active


def process_vcf(vcf_in, vcf_out, ps_active):
    with gzip.open(vcf_in, "rt") as fin, open(vcf_out, "w") as fout:

        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            parts = line.strip().split("\t")

            fmt = parts[8].split(":")
            sample = parts[9].split(":")

            if "PS" not in fmt or "GT" not in fmt:
                fout.write(line)
                continue

            ps = sample[fmt.index("PS")]
            gt = sample[fmt.index("GT")]

            # Skip missing
            if ps not in ps_active or gt == "./.":
                fout.write(line)
                continue

            active = ps_active[ps]

            # Convert to phased genotype
            if gt in ["0/1", "1/0"]:
                if active == "H1":
                    new_gt = "0|1"
                else:
                    new_gt = "1|0"
            else:
                new_gt = gt  # homozygous stays same

            sample[fmt.index("GT")] = new_gt
            parts[9] = ":".join(sample)

            fout.write("\t".join(parts) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Skew-aware phasing")
    parser.add_argument("--vcf", required=True, help="Input phased VCF (gz)")
    parser.add_argument("--skew", required=True, help="Skew TSV (gz)")
    parser.add_argument("--out", required=True, help="Output VCF")

    args = parser.parse_args()

    ps_active = load_skew(args.skew)
    process_vcf(args.vcf, args.out, ps_active)


if __name__ == "__main__":
    main()