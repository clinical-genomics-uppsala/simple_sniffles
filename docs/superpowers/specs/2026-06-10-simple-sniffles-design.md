# simple_sniffles — Design Spec

**Date:** 2026-06-10
**Author:** Andrei Guliaev

## Purpose

A standalone Snakemake pipeline, built to hydra-genetics standards, that takes
already-mapped BAM files (listed in `units.tsv` / `samples.tsv`), runs
**Sniffles2** for structural-variant calling, and produces a per-sample Excel
(`.xlsx`) spreadsheet of the called variants.

Scope is deliberately small: no alignment, no SNV/CNV calling. Input BAMs are
assumed mapped and indexed.

## Inputs

### `config/samples.tsv`
| column | required | notes |
|--------|----------|-------|
| sample | yes | sample id |

### `config/units.tsv`
| column | required | notes |
|--------|----------|-------|
| sample   | yes | matches samples.tsv |
| type     | yes | `N`, `T`, or `R` |
| platform | yes | e.g. `PACBIO`, `ONT` |
| bam      | yes | path to mapped BAM |
| bai      | no  | path to BAM index; defaults to `<bam>.bai` if absent |

Note: alignment-oriented unit keys (`flowcell`, `lane`, `barcode`) are dropped
since the pipeline does not align. The units index becomes `(sample, type)`.

### `config/config.yaml` (key entries)
- `samples`, `units`, `resources`, `output` — paths to the above files.
- `default_container`, `default_resources`.
- `reference.fasta` — reference genome FASTA.
- `reference.tandem_repeats` — tandem-repeat BED passed to Sniffles2.
- `module_versions.cnv_sv` — pinned hydra-genetics/cnv_sv tag.
- Per-rule blocks: `sniffles2_call`, `compile_xlsx_report` (container, resources,
  and `compile_xlsx_report.filters` / `.genes`).

## Architecture

```
config/                       workflow/
  config.yaml                   Snakefile          # rule all + cnv_sv import + use rules
  resources.yaml                rules/common.smk   # load+validate, get_bam() helper
  samples.tsv                   scripts/compile_xlsx_report_sv.py
  units.tsv                     schemas/*.schema.yaml
  output_files.yaml
results/                      # final copied outputs
```

### Data flow
1. **common.smk** loads and validates `samples.tsv` and `units.tsv` against
   schemas; reads `output_files.yaml`. Provides `get_bam(wildcards)` returning
   `(bam, bai)` from the units table (bai derived as `<bam>.bai` when the column
   is empty).
2. **`cnv_sv_sniffles2`** — `use rule sniffles2_call from cnv_sv` (github import,
   tag from `module_versions.cnv_sv`). Inputs: `bam`/`bai` from `get_bam`,
   `ref=reference.fasta`, and `--tandem-repeats reference.tandem_repeats`.
   Output: `cnv_sv/sniffles2_call/{sample}_{type}.vcf`.
3. **`cnv_sv_bgzip_sniffles` / `cnv_sv_tabix_sniffles`** — `use rule bgzip` and
   `use rule tabix from cnv_sv`, constrained to `^cnv_sv/sniffles2_call/.+`,
   producing `{sample}_{type}.vcf.gz` + `.tbi`.
4. **`compile_xlsx_report`** — local rule running
   `scripts/compile_xlsx_report_sv.py`. Input: the bgzipped Sniffles2 VCF.
   Output: `reports/xlsx_reports/{sample}_{type}.sniffles2.xlsx`.
5. **`compile_output_file_list`** (template common.smk) + `output_files.yaml`
   copy the xlsx and the bgzipped VCF into `results/`.

### The xlsx script
`compile_xlsx_report_sv.py` is a trimmed copy of the existing
`examples/workflow/scripts/compile_xlsx_report.py`:
- Keeps SV VCF parsing and `ALLOWED_SV_TYPES = [DEL, INS, INV, DUP, BND]`.
- Drops `vcf_snv` and `vcf_cnv` inputs and their sheets — SV-only.
- `create_report(snakemake_obj)` reads `snakemake.input.vcf_sv`,
  `snakemake.output.xlsx`, `snakemake.params` (`filter_config`, `genes_bed`,
  `software_versions`), and `snakemake.wildcards.sample`.
- Retains the `argparse` fallback for standalone runs (SV args only).

## Rule conventions
Every local rule follows the hydra-genetics rule-template: `log:`, `benchmark:`
with `benchmark_repeats`, `threads`/`resources` resolved via
`config.get(<rule>, {}).get(..., config["default_resources"][...])`,
`container:` with `default_container` fallback, and a `message:`.

## Wildcards
`wildcard_constraints: sample="|".join(samples.index)`, `type="N|T|R"`.

## Schemas
`config`, `resources`, `samples`, `units`, `output_files` schemas adapted from
the pipeline-template. The `units` schema is trimmed to
`sample, type, platform, bam, bai` with `bam` required.

## Error handling
- Config/sample/unit validation via `snakemake.utils.validate` against schemas
  (fails fast with a trimmed message, as in template common.smk).
- Missing BAM path → schema/`get_bam` error before any job runs.
- Sniffles2 / report failures are captured per-rule in `log:` files.

## Testing
- `.tests/integration/` with tiny `samples.tsv`/`units.tsv` pointing at a small
  BAM, plus a `config.yaml`, exercised with `snakemake -n` (dry run) and a full
  run on the test data in CI.

## Out of scope
Alignment, SNV/CNV calling, annotation, multi-sample joint calling, combined
workbooks.
