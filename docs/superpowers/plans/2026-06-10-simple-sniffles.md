# simple_sniffles Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a standalone hydra-genetics Snakemake pipeline that runs Sniffles2 on mapped BAMs listed in `units.tsv`/`samples.tsv` and emits a per-sample xlsx of structural variants.

**Architecture:** A thin workflow that imports the `hydra-genetics/cnv_sv` module via `github()` and reuses its `sniffles2_call`, `bgzip`, and `tabix` rules. A local `common.smk` loads/validates the TSVs and resolves each BAM path from a `bam` column. A local rule runs an SV-only fork of the existing `compile_xlsx_report.py`. Final files are copied to `results/` via the template's `output_files.yaml` machinery.

**Tech Stack:** Snakemake, hydra-genetics utils, pandas, pysam, openpyxl, pytest.

---

## File Structure

- Create: `config/config.yaml` — pipeline config (paths, references, containers, module versions).
- Create: `config/resources.yaml` — default resources.
- Create: `config/samples.tsv` — example samples table.
- Create: `config/units.tsv` — example units table with `bam` column.
- Create: `config/output_files.yaml` — final output spec.
- Create: `workflow/Snakefile` — rule all + cnv_sv module import + `use rule`s + local report rule.
- Create: `workflow/rules/common.smk` — load/validate, `get_bam()` helper, copy-rule machinery.
- Create: `workflow/scripts/compile_xlsx_report_sv.py` — SV-only xlsx generator.
- Create: `workflow/schemas/config.schema.yaml`, `resources.schema.yaml`, `samples.schema.yaml`, `units.schema.yaml`, `output_files.schema.yaml`.
- Create: `.tests/integration/samples.tsv`, `.tests/integration/units.tsv`, `.tests/integration/config.yaml`.
- Create: `tests/test_compile_xlsx_report_sv.py` — unit tests for the xlsx script's parser.
- Create: `requirements.txt`, `README.md`.

---

## Task 1: Scaffold repo and Python requirements

**Files:**
- Create: `requirements.txt`
- Create: `README.md`

- [ ] **Step 1: Write `requirements.txt`**

```
snakemake>=7.32.0,<8
hydra-genetics>=3.0.0
pandas>=2.0
pysam>=0.22
openpyxl>=3.1
PyYAML>=6.0
jsonschema>=4.0
```

- [ ] **Step 2: Write `README.md`**

```markdown
# simple_sniffles

A hydra-genetics Snakemake pipeline that runs Sniffles2 on mapped BAM files
listed in `units.tsv`/`samples.tsv` and produces a per-sample Excel report of
structural variants.

## Run

    snakemake --configfile config/config.yaml --use-singularity -j1

Inputs are mapped BAMs referenced by the `bam` column in `config/units.tsv`.
```

- [ ] **Step 3: Commit**

```bash
git add requirements.txt README.md
git commit -m "chore: scaffold simple_sniffles repo"
```

---

## Task 2: Schemas

**Files:**
- Create: `workflow/schemas/config.schema.yaml`
- Create: `workflow/schemas/resources.schema.yaml`
- Create: `workflow/schemas/samples.schema.yaml`
- Create: `workflow/schemas/units.schema.yaml`
- Create: `workflow/schemas/output_files.schema.yaml`

- [ ] **Step 1: Write `workflow/schemas/config.schema.yaml`**

```yaml
$schema: "http://json-schema.org/draft-04/schema#"
description: snakemake configuration file
type: object
properties:
  samples:
    type: string
  units:
    type: string
  output:
    type: string
    description: output yaml/json file defining expected output from pipeline
  resources:
    type: string
    description: Path to resources.yaml file
  default_container:
    type: string
    description: name or path to a default docker/singularity container
  reference:
    type: object
    properties:
      fasta:
        type: string
        description: path to reference genome FASTA
    required:
      - fasta
  sniffles2_call:
    type: object
    properties:
      container:
        type: string
      tandem_repeats:
        type: string
        description: path to tandem-repeat BED passed to Sniffles2
required:
  - samples
  - units
  - output
  - resources
  - default_container
  - reference
```

- [ ] **Step 2: Write `workflow/schemas/resources.schema.yaml`**

```yaml
$schema: "http://json-schema.org/draft-04/schema#"
description: resource definition
type: object
properties:
  default_resources:
    type: object
    properties:
      mem_mb:
        type: integer
      mem_per_cpu:
        type: integer
      partition:
        type: string
      threads:
        type: integer
      time:
        type: string
required:
  - default_resources
```

- [ ] **Step 3: Write `workflow/schemas/samples.schema.yaml`**

```yaml
$schema: "http://json-schema.org/draft-04/schema#"
description: sample list
properties:
  sample:
    type: string
    description: sample id
required:
  - sample
```

- [ ] **Step 4: Write `workflow/schemas/units.schema.yaml`**

```yaml
$schema: "http://json-schema.org/draft-07/schema#"
description: row represents one mapped BAM dataset
properties:
  sample:
    type: string
    description: sample id
  type:
    type: string
    description: type of sample data Tumor, Normal, RNA (N|T|R)
    pattern: "^(N|T|R)$"
  platform:
    type: string
    description: sequence platform used to generate data
  bam:
    type: string
    description: absolute path to mapped bam file
  bai:
    type: string
    description: absolute path to bam index (optional; defaults to <bam>.bai)
required:
  - sample
  - type
  - platform
  - bam
```

- [ ] **Step 5: Write `workflow/schemas/output_files.schema.yaml`**

```yaml
$schema: "http://json-schema.org/draft-04/schema#"
description: Output file specification
type: object
properties:
  files:
    description: Defines a list of output files
    type: array
    items:
      type: object
      properties:
        name:
          type: string
        input:
          type:
            - string
            - "null"
        output:
          type: string
        types:
          type: array
          items:
            type: string
            pattern: ^(N|T|R)$
      required:
        - name
        - input
        - output
required:
  - files
```

- [ ] **Step 6: Commit**

```bash
git add workflow/schemas/
git commit -m "feat: add config/resources/samples/units/output schemas"
```

---

## Task 3: Config and example TSVs

**Files:**
- Create: `config/config.yaml`
- Create: `config/resources.yaml`
- Create: `config/samples.tsv`
- Create: `config/units.tsv`
- Create: `config/output_files.yaml`

- [ ] **Step 1: Write `config/config.yaml`**

```yaml
---
output: "config/output_files.yaml"
resources: "config/resources.yaml"
samples: "config/samples.tsv"
units: "config/units.tsv"

default_container: "docker://hydragenetics/common:3.0.0"

module_versions:
  cnv_sv: "v2.0.0"

reference:
  fasta: "reference/genome.fasta"

sniffles2_call:
  container: "docker://hydragenetics/sniffles2:2.2"
  tandem_repeats: "reference/human_GRCh38_no_alt_analysis_set.trf.bed"

compile_xlsx_report:
  container: "docker://hydragenetics/common:3.0.0"
  filters: []
  genes: ""
```

- [ ] **Step 2: Write `config/resources.yaml`**

```yaml
---
default_resources:
  threads: 1
  time: "12:00:00"
  mem_mb: 1024
  mem_per_cpu: 3072
  partition: core

sniffles2_call:
  threads: 4
  mem_mb: 8192
```

- [ ] **Step 3: Write `config/samples.tsv`** (tab-separated)

```
sample
HG002
```

- [ ] **Step 4: Write `config/units.tsv`** (tab-separated)

```
sample	type	platform	bam	bai
HG002	T	PACBIO	/data/HG002_T.bam	/data/HG002_T.bam.bai
```

- [ ] **Step 5: Write `config/output_files.yaml`**

```yaml
---
directory: "results"
files:
  - name: "sniffles2 vcf"
    input: "cnv_sv/sniffles2_call/{sample}_{type}.vcf.gz"
    output: "sv/{sample}_{type}.sniffles2.vcf.gz"
  - name: "sniffles2 xlsx report"
    input: "reports/xlsx_reports/{sample}_{type}.sniffles2.xlsx"
    output: "sv/{sample}_{type}.sniffles2.xlsx"
```

- [ ] **Step 6: Commit**

```bash
git add config/
git commit -m "feat: add config, resources, example TSVs and output spec"
```

---

## Task 4: common.smk

**Files:**
- Create: `workflow/rules/common.smk`

- [ ] **Step 1: Write `workflow/rules/common.smk`**

Adapt the pipeline-template `common.smk`. Key changes vs template: units index is `(sample, type)`; add a `get_bam(wildcards)` helper; keep `compile_output_file_list` and `generate_copy_rules` verbatim from the template.

```python
__author__ = "Andrei Guliaev"
__copyright__ = "Copyright 2026, Andrei Guliaev"
__email__ = "andrei.guliaev@scilifelab.uu.se"
__license__ = "GPL-3"

import json
import pathlib
import re
import sys
from pathlib import Path

import pandas as pd
import pandas
import yaml
from snakemake.utils import validate
from snakemake.utils import min_version
from snakemake.exceptions import WorkflowError

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import get_samples
from hydra_genetics.utils.units import get_unit_types

min_version("7.32.0")

### Validate config

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")

validate(config, schema="../schemas/config.schema.yaml")

### Resources

config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")

### Samples

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Units (indexed by sample + type only)

units = (
    pandas.read_table(config["units"], dtype=str)
    .set_index(["sample", "type"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")

### Output spec

with open(config["output"]) as f:
    if config["output"].endswith("json"):
        output_spec = json.load(f)
    else:
        output_spec = yaml.safe_load(f.read())
validate(output_spec, schema="../schemas/output_files.schema.yaml")

### Wildcards

wildcard_constraints:
    sample="|".join(samples.index),
    type="N|T|R",


def get_bam(wildcards):
    """Return (bam, bai) for a sample/type from units.tsv.

    bai falls back to <bam>.bai when the column is empty/missing.
    """
    unit = units.loc[(wildcards.sample, wildcards.type)]
    bam = unit["bam"]
    bai = unit.get("bai", "")
    if not isinstance(bai, str) or bai == "" or pd.isna(bai):
        bai = f"{bam}.bai"
    return bam, bai


def compile_output_file_list(wildcards):
    outdir = pathlib.Path(output_spec.get("directory", "./"))
    output_files = []
    for f in output_spec["files"]:
        outputpaths = set(
            [
                f["output"].format(sample=sample, type=unit_type)
                for sample in get_samples(samples)
                for unit_type in get_unit_types(units, sample)
            ]
        )
        for op in outputpaths:
            output_files.append(outdir / Path(op))
    return output_files


def generate_copy_rules(output_spec):
    output_directory = pathlib.Path(output_spec.get("directory", "./"))
    rulestrings = []
    for f in output_spec["files"]:
        if f["input"] is None:
            continue
        rule_name = "_copy_{}".format("_".join(re.split(r"\W{1,}", f["name"].strip().lower())))
        input_file = pathlib.Path(f["input"])
        output_file = output_directory / pathlib.Path(f["output"])
        mem_mb = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
        mem_per_cpu = config.get("_copy", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"])
        partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
        threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
        time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
        copy_container = config.get("_copy", {}).get("container", config["default_container"])
        rule_code = "\n".join(
            [
                f'@workflow.rule(name="{rule_name}")',
                f'@workflow.input("{input_file}")',
                f'@workflow.output("{output_file}")',
                f'@workflow.log("logs/{rule_name}_{output_file.name}.log")',
                f'@workflow.container("{copy_container}")',
                f'@workflow.resources(time="{time}", threads={threads}, mem_mb="{mem_mb}", '
                f'mem_per_cpu={mem_per_cpu}, partition="{partition}")',
                f'@workflow.shellcmd("{copy_container}")',
                "@workflow.run\n",
                f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, "
                "log, version, rule, conda_env, container_img, singularity_args, use_singularity, "
                "env_modules, bench_record, jobid, is_shell, bench_iteration, cleanup_scripts, "
                "shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
                "__is_snakemake_rule_func=True):",
                '\tshell("(cp --preserve=timestamps {input[0]} {output[0]}) &> {log}", bench_record=bench_record, '
                "bench_iteration=bench_iteration)\n\n",
            ]
        )
        rulestrings.append(rule_code)
    exec(compile("\n".join(rulestrings), "copy_result_files", "exec"), workflow.globals)


generate_copy_rules(output_spec)
```

- [ ] **Step 2: Commit**

```bash
git add workflow/rules/common.smk
git commit -m "feat: add common.smk with get_bam helper and copy rules"
```

---

## Task 5: SV-only xlsx report script — failing test first

**Files:**
- Create: `tests/test_compile_xlsx_report_sv.py`
- Create: `workflow/scripts/compile_xlsx_report_sv.py` (stub in this task)

- [ ] **Step 1: Create a tiny Sniffles2-style VCF fixture and write the failing test**

Create `tests/test_compile_xlsx_report_sv.py`:

```python
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
```

- [ ] **Step 2: Create a stub script so the import resolves**

Create `workflow/scripts/compile_xlsx_report_sv.py`:

```python
def parse_sv_vcf(path):
    raise NotImplementedError


def write_xlsx(rows, out_path, sample, software_versions=None):
    raise NotImplementedError
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `pytest tests/test_compile_xlsx_report_sv.py -v`
Expected: FAIL with `NotImplementedError`.

- [ ] **Step 4: Commit**

```bash
git add tests/test_compile_xlsx_report_sv.py workflow/scripts/compile_xlsx_report_sv.py
git commit -m "test: failing tests for SV-only xlsx report"
```

---

## Task 6: Implement the SV-only xlsx report script

**Files:**
- Modify: `workflow/scripts/compile_xlsx_report_sv.py`

- [ ] **Step 1: Implement the full script**

Replace `workflow/scripts/compile_xlsx_report_sv.py` with:

```python
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
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
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
```

- [ ] **Step 2: Run the tests to verify they pass**

Run: `pytest tests/test_compile_xlsx_report_sv.py -v`
Expected: PASS (2 passed).

- [ ] **Step 3: Commit**

```bash
git add workflow/scripts/compile_xlsx_report_sv.py
git commit -m "feat: implement SV-only xlsx report script"
```

---

## Task 7: Snakefile — module import, use rules, report rule

**Files:**
- Create: `workflow/Snakefile`

- [ ] **Step 1: Write `workflow/Snakefile`**

```python
__author__ = "Andrei Guliaev"
__copyright__ = "Copyright 2026, Andrei Guliaev"
__email__ = "andrei.guliaev@scilifelab.uu.se"
__license__ = "GPL-3"


include: "rules/common.smk"


rule all:
    input:
        compile_output_file_list,


module cnv_sv:
    snakefile:
        github(
            "hydra-genetics/cnv_sv",
            path="workflow/Snakefile",
            tag=config.get("module_versions", {}).get("cnv_sv", "v2.0.0"),
        )
    config:
        config


# Run Sniffles2 on the mapped BAM referenced in units.tsv.
# Override bam/bai/ref; tandem-repeats comes from config["sniffles2_call"]["tandem_repeats"]
# via the module's get_tr_bed (params.tandem_repeats), so no params override is needed.
use rule sniffles2_call from cnv_sv as cnv_sv_sniffles2 with:
    input:
        bam=lambda wildcards: get_bam(wildcards)[0],
        bai=lambda wildcards: get_bam(wildcards)[1],
        ref=config.get("reference", {}).get("fasta", ""),


# bgzip + tabix the Sniffles2 VCF
use rule bgzip from cnv_sv as cnv_sv_bgzip_sniffles with:
    wildcard_constraints:
        file="^cnv_sv/sniffles2_call/.+",


use rule tabix from cnv_sv as cnv_sv_tabix_sniffles with:
    wildcard_constraints:
        file="^cnv_sv/sniffles2_call/.+",


rule compile_xlsx_report:
    input:
        vcf_sv="cnv_sv/sniffles2_call/{sample}_{type}.vcf.gz",
        tbi="cnv_sv/sniffles2_call/{sample}_{type}.vcf.gz.tbi",
    output:
        xlsx=temp("reports/xlsx_reports/{sample}_{type}.sniffles2.xlsx"),
    params:
        software_versions={
            "sniffles2": config.get("sniffles2_call", {}).get("container", ""),
        },
    log:
        "reports/xlsx_reports/{sample}_{type}.sniffles2.xlsx.log",
    benchmark:
        repeat(
            "reports/xlsx_reports/{sample}_{type}.sniffles2.xlsx.benchmark.tsv",
            config.get("compile_xlsx_report", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("compile_xlsx_report", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("compile_xlsx_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("compile_xlsx_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("compile_xlsx_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("compile_xlsx_report", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("compile_xlsx_report", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("compile_xlsx_report", {}).get("container", config["default_container"])
    message:
        "{rule}: compile Sniffles2 SVs into Excel report for {wildcards.sample}_{wildcards.type}"
    script:
        "scripts/compile_xlsx_report_sv.py"
```

- [ ] **Step 2: Commit**

```bash
git add workflow/Snakefile
git commit -m "feat: add Snakefile with cnv_sv import and report rule"
```

---

## Task 8: Integration test config and dry-run validation

**Files:**
- Create: `.tests/integration/samples.tsv`
- Create: `.tests/integration/units.tsv`
- Create: `.tests/integration/config.yaml`

- [ ] **Step 1: Write `.tests/integration/samples.tsv`** (tab-separated)

```
sample
testsample
```

- [ ] **Step 2: Write `.tests/integration/units.tsv`** (tab-separated)

```
sample	type	platform	bam	bai
testsample	T	PACBIO	.tests/integration/data/testsample_T.bam	.tests/integration/data/testsample_T.bam.bai
```

- [ ] **Step 3: Write `.tests/integration/config.yaml`**

```yaml
---
output: "config/output_files.yaml"
resources: "config/resources.yaml"
samples: ".tests/integration/samples.tsv"
units: ".tests/integration/units.tsv"

default_container: "docker://hydragenetics/common:3.0.0"

module_versions:
  cnv_sv: "v2.0.0"

reference:
  fasta: ".tests/integration/data/reference.fasta"

sniffles2_call:
  container: "docker://hydragenetics/sniffles2:2.2"
  tandem_repeats: ".tests/integration/data/tandem_repeats.bed"

compile_xlsx_report:
  container: "docker://hydragenetics/common:3.0.0"
```

- [ ] **Step 4: Dry-run the workflow to validate wiring**

Run: `snakemake -n --configfile .tests/integration/config.yaml -s workflow/Snakefile`
Expected: Snakemake builds a DAG with `cnv_sv_sniffles2`, `cnv_sv_bgzip_sniffles`, `cnv_sv_tabix_sniffles`, `compile_xlsx_report`, and `_copy_*` jobs for `testsample_T`, and prints "Job stats" without config/schema validation errors. (Requires network access to fetch the cnv_sv module from GitHub.)

- [ ] **Step 5: Commit**

```bash
git add .tests/integration/
git commit -m "test: add integration config and dry-run target"
```

---

## Task 9: Final verification

- [ ] **Step 1: Run unit tests**

Run: `pytest tests/ -v`
Expected: all pass.

- [ ] **Step 2: Lint the Snakefile syntax via dry run**

Run: `snakemake -n --configfile .tests/integration/config.yaml -s workflow/Snakefile`
Expected: DAG resolves, no errors.

- [ ] **Step 3: Commit any fixups**

```bash
git add -A
git commit -m "chore: final verification fixups" || echo "nothing to commit"
```

---

## Notes for the implementer

- The `cnv_sv` `sniffles2_call` rule reads tandem repeats from
  `config["sniffles2_call"]["tandem_repeats"]` via its own `get_tr_bed`; do NOT
  redefine `get_tr_bed` locally — setting the config key is enough.
- The module's default `bam` input uses `get_input_aligned_bam`; we fully
  override `bam`/`bai` with `get_bam`, so that helper is never invoked.
- `sniffles2_call` also emits a `.snf` output (kept as the module defines it);
  we only consume the `.vcf`.
- Container tags (`hydragenetics/sniffles2:2.2`, `common:3.0.0`) and the
  `cnv_sv` module tag (`v2.0.0`) are reasonable defaults — adjust to match your
  deployed versions before a production run.
```
