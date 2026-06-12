"""Compile a Sniffles2 structural-variant VCF into an xlsx report.

Runs either as a Snakemake script (uses the injected `snakemake` object) or
standalone via argparse. Handles both single-sample and joint (multi-sample)
Sniffles2 VCFs.
"""
import argparse
import json
import logging
from typing import Any

import pandas as pd
import pysam

ALLOWED_SV_TYPES = ["DEL", "INS", "INV", "DUP", "BND"]

# Column order for single-sample output; joint output extends this dynamically.
SINGLE_SAMPLE_COLUMNS = [
    "CHROM", "POS", "ID", "SVTYPE", "END", "SVLEN",
    "ALT", "QUAL", "FILTER", "GT", "DR", "DV", "VAF", "COVERAGE",
]


def parse_sv_vcf(path: str, sample_names: list = None) -> list:
    """Parse a Sniffles2 VCF with pysam, keeping allowed SV types.

    Args:
        path: Path to VCF (plain or bgzipped).
        sample_names: Sample column names to extract FORMAT fields from.
            None (default) uses the first sample in the header → single-sample
            mode with flat column names (GT, DR, DV, VAF).
            A list of names activates multi-sample mode; FORMAT columns are
            suffixed with the type letter from the sample name, e.g. GT_T, GT_N.
    """
    rows = []
    with pysam.VariantFile(path, "r") as vcf:
        all_samples = list(vcf.header.samples)
        if sample_names is None:
            sample_names = all_samples[:1]

        # In multi-sample mode suffix columns with the type letter (last _-component).
        multi = len(sample_names) > 1
        suffixes = {s: f"_{s.rsplit('_', 1)[-1]}" if multi else "" for s in sample_names}

        for record in vcf:
            svtype = record.info["SVTYPE"] if "SVTYPE" in record.info else ""
            if svtype not in ALLOWED_SV_TYPES:
                continue

            svlen = record.info["SVLEN"] if "SVLEN" in record.info else None
            if isinstance(svlen, (list, tuple)):
                svlen = svlen[0] if svlen else None

            raw_cov = record.info["COVERAGE"] if "COVERAGE" in record.info else None
            if raw_cov is None:
                coverage = ""
            elif isinstance(raw_cov, (list, tuple)):
                coverage = ",".join(str(int(v)) for v in raw_cov if v is not None)
            else:
                coverage = str(raw_cov)

            row = {
                "CHROM": record.chrom,
                "POS": record.pos,
                "ID": record.id or "",
                "SVTYPE": svtype,
                "END": record.stop,
                "SVLEN": svlen,
                "ALT": record.alts[0] if record.alts else "",
                "QUAL": record.qual,
                "FILTER": ";".join(record.filter.keys()) or ".",
            }

            for sname in sample_names:
                sfx = suffixes[sname]
                samp = record.samples[sname] if sname in record.samples else None

                gt = ""
                if samp is not None:
                    raw_gt = samp.get("GT")
                    if raw_gt is not None:
                        gt = "/".join("." if a is None else str(a) for a in raw_gt)

                dr = samp.get("DR") if samp is not None else None
                dv = samp.get("DV") if samp is not None else None

                # VAF: INFO first (Sniffles2 writes it there), FORMAT fallback
                vaf = record.info["VAF"] if "VAF" in record.info else None
                if vaf is None and samp is not None:
                    vaf = samp.get("VAF")
                if isinstance(vaf, (list, tuple)):
                    vaf = vaf[0] if vaf else None

                row[f"GT{sfx}"] = gt
                row[f"DR{sfx}"] = dr
                row[f"DV{sfx}"] = dv
                row[f"VAF{sfx}"] = vaf

            row["COVERAGE"] = coverage
            rows.append(row)
    return rows


def write_xlsx(rows: list, out_path: str, sample: str, software_versions: dict = None) -> None:
    software_versions = software_versions or {}
    columns = list(rows[0].keys()) if rows else SINGLE_SAMPLE_COLUMNS
    variants = pd.DataFrame(rows, columns=columns)
    meta = pd.DataFrame(
        [{"key": "sample", "value": sample}]
        + [{"key": f"version:{k}", "value": v} for k, v in software_versions.items()]
    )
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        variants.to_excel(writer, sheet_name="SV", index=False)
        meta.to_excel(writer, sheet_name="info", index=False)


def create_report(snakemake_obj: Any) -> None:
    log_file = snakemake_obj.log[0] if snakemake_obj.log else None
    logging.basicConfig(filename=log_file, level=logging.INFO)
    vcf_sv = snakemake_obj.input.vcf_sv
    out_xlsx = snakemake_obj.output.xlsx
    sample = snakemake_obj.wildcards.sample
    sample_names = getattr(snakemake_obj.params, "sample_names", None)
    software_versions = getattr(snakemake_obj.params, "software_versions", {})
    logging.info(f"Parsing {vcf_sv} for sample {sample} (sample_names={sample_names})")
    rows = parse_sv_vcf(vcf_sv, sample_names=sample_names)
    write_xlsx(rows, out_xlsx, sample=sample, software_versions=software_versions)
    logging.info(f"Wrote {len(rows)} SV rows to {out_xlsx}")


if __name__ == "__main__":
    if "snakemake" in dir():
        create_report(snakemake)  # type: ignore  # noqa: F821
    else:
        parser = argparse.ArgumentParser(description="Compile SV xlsx report from a Sniffles2 VCF.")
        parser.add_argument("--vcf-sv", required=True)
        parser.add_argument("--output-xlsx", "-o", required=True)
        parser.add_argument("--sample", default="unknown")
        parser.add_argument("--sample-names", nargs="+", default=None,
                            help="Sample column names to extract (multi-sample / joint mode).")
        parser.add_argument("--software-versions", default="{}")
        args = parser.parse_args()
        write_xlsx(
            parse_sv_vcf(args.vcf_sv, sample_names=args.sample_names),
            args.output_xlsx,
            sample=args.sample,
            software_versions=json.loads(args.software_versions),
        )
