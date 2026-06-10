"""Compile a Sniffles2 structural-variant VCF into an xlsx report.

Runs either as a Snakemake script (uses the injected `snakemake` object) or
standalone via argparse. SV-only: parses one Sniffles2 VCF and writes one
spreadsheet.
"""
import argparse
import json
import logging
from typing import Any

import pandas as pd
import pysam

ALLOWED_SV_TYPES = ["DEL", "INS", "INV", "DUP", "BND"]

COLUMNS = [
    "CHROM", "POS", "ID", "SVTYPE", "END", "SVLEN",
    "ALT", "QUAL", "FILTER", "GT", "DR", "DV", "VAF", "COVERAGE",
]


def parse_sv_vcf(path: str) -> list:
    """Parse a Sniffles2 VCF with pysam, keeping allowed SV types.

    Extracts FORMAT fields GT, DR, DV, VAF and INFO field COVERAGE.
    """
    rows = []
    with pysam.VariantFile(path, "r") as vcf:
        samples = list(vcf.header.samples)
        samp_name = samples[0] if samples else None

        for record in vcf:
            svtype = record.info["SVTYPE"] if "SVTYPE" in record.info else ""
            if svtype not in ALLOWED_SV_TYPES:
                continue

            samp = record.samples[samp_name] if samp_name else None

            # Genotype
            gt = ""
            if samp is not None:
                raw_gt = samp.get("GT")
                if raw_gt is not None:
                    gt = "/".join("." if a is None else str(a) for a in raw_gt)

            # Read support
            dr = samp.get("DR") if samp is not None else None
            dv = samp.get("DV") if samp is not None else None

            # VAF: INFO first (Sniffles2 writes it there), FORMAT fallback
            vaf = record.info["VAF"] if "VAF" in record.info else None
            if vaf is None and samp is not None:
                vaf = samp.get("VAF")
            if isinstance(vaf, (list, tuple)):
                vaf = vaf[0] if vaf else None

            # COVERAGE: Sniffles2 emits a tuple; join as comma-separated string
            raw_cov = record.info["COVERAGE"] if "COVERAGE" in record.info else None
            if raw_cov is None:
                coverage = ""
            elif isinstance(raw_cov, (list, tuple)):
                coverage = ",".join(str(int(v)) for v in raw_cov if v is not None)
            else:
                coverage = str(raw_cov)

            svlen = record.info["SVLEN"] if "SVLEN" in record.info else None
            if isinstance(svlen, (list, tuple)):
                svlen = svlen[0] if svlen else None

            rows.append(
                {
                    "CHROM": record.chrom,
                    "POS": record.pos,
                    "ID": record.id or "",
                    "SVTYPE": svtype,
                    "END": record.stop,
                    "SVLEN": svlen,
                    "ALT": record.alts[0] if record.alts else "",
                    "QUAL": record.qual,
                    "FILTER": ";".join(record.filter.keys()) or ".",
                    "GT": gt,
                    "DR": dr,
                    "DV": dv,
                    "VAF": vaf,
                    "COVERAGE": coverage,
                }
            )
    return rows


def write_xlsx(rows: list, out_path: str, sample: str, software_versions: dict = None) -> None:
    software_versions = software_versions or {}
    variants = pd.DataFrame(rows, columns=COLUMNS)
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
