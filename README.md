![CodeRabbit Pull Request Reviews](https://img.shields.io/coderabbit/prs/github/clinical-genomics-uppsala/simple_sniffles?utm_source=oss&utm_medium=github&utm_campaign=clinical-genomics-uppsala%2Fsimple_sniffles&labelColor=171717&color=FF570A&link=https%3A%2F%2Fcoderabbit.ai&label=CodeRabbit+Reviews)

# A simple Sniffles2 pipeline

A hydra-genetics Snakemake pipeline that runs Sniffles2 on mapped BAM files
listed in `units.tsv`/`samples.tsv` and produces a per-sample Excel report of
structural variants.

## Setup

    pixi install

## Run

    pixi run snakemake --configfile config/config.yaml --use-singularity -j1

Inputs are mapped BAMs referenced by the `bam` column in `config/units.tsv` - the file was created manually not with help of `hydra-genetics create-input-files` because they are mapped files.
