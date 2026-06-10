import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).parent.parent / "workflow" / "scripts" / "compile_xlsx_report_sv.py"

VCF = """\
##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002
chr1\t1000\tSniffles2.DEL.1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=-500;END=1500\tGT\t0/1
chr2\t2000\tSniffles2.INS.1\tN\t<INS>\t60\tPASS\tSVTYPE=INS;SVLEN=300;END=2000\tGT\t1/1
chr3\t3000\tSniffles2.XYZ.1\tN\t<XYZ>\t60\tPASS\tSVTYPE=XYZ;END=3100\tGT\t0/1
"""


def load_module():
    spec = importlib.util.spec_from_file_location("xlsx_sv", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_parse_sv_vcf_keeps_allowed_types(tmp_path):
    vcf_path = tmp_path / "sv.vcf"
    vcf_path.write_text(VCF)
    mod = load_module()
    rows = mod.parse_sv_vcf(str(vcf_path))
    svtypes = sorted(r["SVTYPE"] for r in rows)
    assert svtypes == ["DEL", "INS"]  # XYZ filtered out
    del_row = next(r for r in rows if r["SVTYPE"] == "DEL")
    assert del_row["CHROM"] == "chr1"
    assert del_row["POS"] == 1000
    assert del_row["SVLEN"] == -500


def test_write_xlsx_creates_file(tmp_path):
    vcf_path = tmp_path / "sv.vcf"
    vcf_path.write_text(VCF)
    out = tmp_path / "out.xlsx"
    mod = load_module()
    rows = mod.parse_sv_vcf(str(vcf_path))
    mod.write_xlsx(rows, str(out), sample="HG002", software_versions={"sniffles2": "2.2"})
    assert out.exists() and out.stat().st_size > 0
