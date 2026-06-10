import importlib.util
from pathlib import Path

import openpyxl

SCRIPT = Path(__file__).parent.parent / "workflow" / "scripts" / "compile_xlsx_report_sv.py"

VCF_LINES = """\
##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">
##INFO=<ID=COVERAGE,Number=.,Type=Integer,Description="Coverage">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Ref reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Alt reads">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002
chr1\t999\tSniffles2.DEL.1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=-500;END=1500;VAF=0.5;COVERAGE=30,28\tGT:DR:DV\t0/1:15:15
chr2\t1999\tSniffles2.INS.1\tN\t<INS>\t60\tPASS\tSVTYPE=INS;SVLEN=300;END=2000;VAF=0.9;COVERAGE=20\tGT:DR:DV\t1/1:2:18
chr3\t2999\tSniffles2.XYZ.1\tN\t<XYZ>\t60\tPASS\tSVTYPE=XYZ;END=3100\tGT:DR:DV\t0/1:10:5
"""


def make_vcf(tmp_path):
    """Write a plain VCF; pysam.VariantFile reads it without bgzip/tabix."""
    vcf = tmp_path / "sv.vcf"
    vcf.write_text(VCF_LINES)
    return str(vcf)


def load_module():
    spec = importlib.util.spec_from_file_location("xlsx_sv", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_parse_sv_vcf_keeps_allowed_types(tmp_path):
    vcf_gz = make_vcf(tmp_path)
    mod = load_module()
    rows = mod.parse_sv_vcf(vcf_gz)
    svtypes = sorted(r["SVTYPE"] for r in rows)
    assert svtypes == ["DEL", "INS"]  # XYZ filtered out


def test_parse_sv_vcf_fields(tmp_path):
    vcf_gz = make_vcf(tmp_path)
    mod = load_module()
    rows = mod.parse_sv_vcf(vcf_gz)
    del_row = next(r for r in rows if r["SVTYPE"] == "DEL")
    assert del_row["CHROM"] == "chr1"
    assert del_row["POS"] == 999      # pysam returns 0-based POS
    assert del_row["SVLEN"] == -500
    assert del_row["GT"] == "0/1"
    assert del_row["DR"] == 15
    assert del_row["DV"] == 15
    assert abs(del_row["VAF"] - 0.5) < 1e-6
    assert del_row["COVERAGE"] == "30,28"


def test_write_xlsx_creates_file(tmp_path):
    vcf_gz = make_vcf(tmp_path)
    out = tmp_path / "out.xlsx"
    mod = load_module()
    rows = mod.parse_sv_vcf(vcf_gz)
    mod.write_xlsx(rows, str(out), sample="HG002", software_versions={"sniffles2": "2.2"})
    assert out.exists() and out.stat().st_size > 0


def test_write_xlsx_columns(tmp_path):
    vcf_gz = make_vcf(tmp_path)
    out = tmp_path / "out.xlsx"
    mod = load_module()
    rows = mod.parse_sv_vcf(vcf_gz)
    mod.write_xlsx(rows, str(out), sample="HG002", software_versions={"sniffles2": "2.2"})
    wb = openpyxl.load_workbook(str(out), read_only=True)
    headers = [cell.value for cell in next(wb["SV"].rows)]
    for col in ("CHROM", "POS", "SVTYPE", "GT", "DR", "DV", "VAF", "COVERAGE"):
        assert col in headers, f"{col} missing from SV sheet"
