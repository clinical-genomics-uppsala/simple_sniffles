# simple_sniffles

A hydra-genetics Snakemake pipeline that runs Sniffles2 on mapped BAM files
listed in `units.tsv`/`samples.tsv` and produces a per-sample Excel report of
structural variants.

## Run

    snakemake --configfile config/config.yaml --use-singularity -j1

Inputs are mapped BAMs referenced by the `bam` column in `config/units.tsv`.
