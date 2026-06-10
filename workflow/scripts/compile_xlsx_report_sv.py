"""Compile a Sniffles2 structural-variant VCF into an xlsx report.

Runs either as a Snakemake script (uses the injected `snakemake` object) or
standalone via argparse. SV-only: parses one Sniffles2 VCF and writes one
spreadsheet.
"""
import argparse
import gzip
import json
import logging
from typing import Any

import pandas as pd

ALLOWED_SV_TYPES = ["DEL", "INS", "INV", "DUP", "BND"]


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _parse_info(info: str) -> dict:
    d = {}
    for field in info.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            d[k] = v
        else:
            d[field] = True
    return d


def parse_sv_vcf(path: str) -> list:
    """Parse a Sniffles2 VCF into a list of row dicts, keeping allowed SV types."""
    rows = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt, qual, flt, info = cols[:8]
            i = _parse_info(info)
            svtype = i.get("SVTYPE", "")
            if svtype not in ALLOWED_SV_TYPES:
                continue
            rows.append(
                {
                    "CHROM": chrom,
                    "POS": int(pos),
                    "ID": vid,
                    "SVTYPE": svtype,
                    "END": int(i["END"]) if "END" in i and i["END"] is not True else None,
                    "SVLEN": int(i["SVLEN"]) if "SVLEN" in i and i["SVLEN"] is not True else None,
                    "ALT": alt,
                    "QUAL": qual,
                    "FILTER": flt,
                }
            )
    return rows


def write_xlsx(rows: list, out_path: str, sample: str, software_versions: dict = None) -> None:
    software_versions = software_versions or {}
    variants = pd.DataFrame(
        rows, columns=["CHROM", "POS", "ID", "SVTYPE", "END", "SVLEN", "ALT", "QUAL", "FILTER"]
    )
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
    software_versions = getattr(snakemake_obj.params, "software_versions", {})
    logging.info(f"Parsing {vcf_sv} for sample {sample}")
    rows = parse_sv_vcf(vcf_sv)
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
        parser.add_argument("--software-versions", default="{}")
        args = parser.parse_args()
        write_xlsx(
            parse_sv_vcf(args.vcf_sv),
            args.output_xlsx,
            sample=args.sample,
            software_versions=json.loads(args.software_versions),
        )
