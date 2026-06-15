# simple_sniffles

A hydra-genetics Snakemake pipeline that runs Sniffles2 on mapped BAM files
listed in `units.tsv`/`samples.tsv` and produces a per-sample Excel report of
structural variants.

## Setup

    pixi install

## Run

    pixi run snakemake --configfile config/config.yaml --use-singularity -j1

Inputs are mapped BAMs referenced by the `bam` column in `config/units.tsv` - the file was created manually not with help of `hydra-genetics create-input-files` because they are mapped files.
