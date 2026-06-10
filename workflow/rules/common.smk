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
