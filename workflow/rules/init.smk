##################################################
# INIT: Common setup for the workflow
##################################################

# ----------------------------
# MODULES
# ----------------------------
import os
import re
import os.path
import pandas
from collections.abc import Mapping
from snakemake.utils import validate

##################################################
# CONFIG VALIDATION + SAMPLES
##################################################

# Validate top-level config
validate(config, srcdir("../../schemas/config.schema.yaml"))

cfg_samples = config.get("samples")

# Branch on type of `samples`: string path vs inline dict with bulk/rhizo
if isinstance(cfg_samples, str):
    # ---- Mode A: `samples` is a path to a TSV/CSV sample table ----
    SAMPLES = pandas.read_csv(cfg_samples, header=0, sep="\t").set_index("Sample_ID", drop=False)
    # Optional: validate the sample table if you keep a schema for it
    try:
        validate(SAMPLES, srcdir("../../schemas/samples.schema.yaml"))
    except Exception:
        # Comment out the try/except if you always keep this schema
        pass

    # Determine available read types from columns
    has_sR1_col = "sR1" in SAMPLES.columns
    has_sR2_col = "sR2" in SAMPLES.columns
    has_long_col = "long_reads" in SAMPLES.columns

    has_short = (
        has_sR1_col and has_sR2_col
        and (SAMPLES["sR1"].notna() & SAMPLES["sR2"].notna()).any()
    )
    has_long = has_long_col and SAMPLES["long_reads"].notna().any()

    READ_TYPES = []
    if has_short:
        READ_TYPES.append("short")
    if has_long:
        READ_TYPES.append("long")

    # If you also want BULK/RHIZO lists from this table, derive here (optional)
    # BULK_SAMPLES  = SAMPLES.loc[SAMPLES["group"] == "bulk",  "Sample_ID"].tolist()
    # RHIZO_SAMPLES = SAMPLES.loc[SAMPLES["group"] == "rhizo", "Sample_ID"].tolist()
    # ALL_HIFI_SAMPLES = BULK_SAMPLES + RHIZO_SAMPLES

    SAMPLES_FILE = os.path.abspath(cfg_samples)

elif isinstance(cfg_samples, Mapping):
    # ---- Mode B: `samples` is an inline dict with bulk/rhizo arrays ----
    BULK_SAMPLES   = list(cfg_samples.get("bulk", []))
    RHIZO_SAMPLES  = list(cfg_samples.get("rhizo", []))
    ALL_HIFI_SAMPLES = BULK_SAMPLES + RHIZO_SAMPLES

    # Fabricate a minimal SAMPLES table for downstream code that expects it
    SAMPLES = pandas.DataFrame({
        "Sample_ID": ALL_HIFI_SAMPLES,
        "group": (["bulk"] * len(BULK_SAMPLES)) + (["rhizo"] * len(RHIZO_SAMPLES)),
    }).set_index("Sample_ID", drop=False)

    # In this inline mode, we assume long-read workflow
    READ_TYPES = ["long"]

    SAMPLES_FILE = None
else:
    raise WorkflowError("config['samples'] must be either a path (str) or a mapping with 'bulk'/'rhizo' arrays.")

# Final safety check
if not READ_TYPES:
    raise WorkflowError(
        "No valid reads detected. If using an inline 'samples' dict, READ_TYPES=['long'] is set automatically; "
        "if using a sample table, ensure sR1/sR2 and/or long_reads columns are present."
    )

##################################################
# PATHS
##################################################

SRC_DIR      = srcdir("../../scripts")        # additional scripts
ENV_DIR      = srcdir("../../envs")           # conda environments
MOD_DIR      = srcdir("../../submodules")     # git submodules

DBS_DIR      = config["dbsdir"]               # database folder (string path as given)
MDATA_DIR    = os.path.abspath(config.get("metadata", "."))  # metadata root (optional)

OLDPWD       = os.path.abspath(os.getcwd())   # original working directory

RESULTS_DIR   = os.path.abspath(config["results_dir"])
LOG_DIR       = os.path.abspath(config["log_dir"])
BENCHMARK_DIR = os.path.abspath(config["benchmark_dir"])

# Sample IDs (only defined here if inline mode was used; if table mode is used
# and you want these lists, derive them above where indicated)
if isinstance(cfg_samples, Mapping):
    # already set earlier
    pass

HYBRID_INDIVIDUAL = config["hybrid_individual_sample"]         # sample ID (not a path)
ONT_FILE          = os.path.abspath(config["ont_file"])         # real path

##################################################
# EXECUTION
##################################################

# default shell for snakemake
shell.executable("bash")

# working directory from config
workdir:
    config["workdir"]

##################################################
# LR HELPERS (HiFi paths)
##################################################

def hifi_path_for_sample(sample: str) -> str:
    """
    Return the HiFi fastq path for a sample.
    Normal: {hifi_dir}/{sample}/{sample}_trimmed.fastq.gz
    Special-case for b26_t3_con:
        {hifi_dir}/{sample}/CEH_hifi_sample_trimmed.fastq.gz
    """
    base = os.path.abspath(config["hifi_dir"])
    p1 = os.path.join(base, sample, f"{sample}_trimmed.fastq.gz")
    p2 = os.path.join(base, sample, "CEH_hifi_sample_trimmed.fastq.gz")
    if os.path.exists(p1):
        return p1
    if os.path.exists(p2):
        return p2
    # Fall back to canonical path so Snakemake raises MissingInputException
    return p1

def get_hifi_reads(wildcards):
    """Single-sample helper used by individual assembly rules."""
    return hifi_path_for_sample(wildcards.sample)

#def get_bulk_hifi_reads(_wildcards=None):
#    """Space-separated HiFi files for bulk co-assembly."""
#    if not isinstance(cfg_samples, Mapping):
#        raise WorkflowError("Bulk co-assembly requested but 'samples' is not an inline dict with 'bulk' entries.")
#    return " ".join(hifi_path_for_sample(s) for s in BULK_SAMPLES)

#def get_rhizo_hifi_reads(_wildcards=None):
#    """Space-separated HiFi files for rhizo co-assembly."""
#    if not isinstance(cfg_samples, Mapping):
#        raise WorkflowError("Rhizo co-assembly requested but 'samples' is not an inline dict with 'rhizo' entries.")
#    return " ".join(hifi_path_for_sample(s) for s in RHIZO_SAMPLES)

def get_bulk_hifi_reads(wildcards):
    """List of HiFi files for bulk co-assembly."""
    if not isinstance(cfg_samples, Mapping):
        raise WorkflowError(
            "Bulk co-assembly requested but 'samples' is not an inline dict with 'bulk' entries."
        )
    return [hifi_path_for_sample(s) for s in BULK_SAMPLES]

def get_rhizo_hifi_reads(wildcards):
    """List of HiFi files for rhizo co-assembly."""
    if not isinstance(cfg_samples, Mapping):
        raise WorkflowError(
            "Rhizo co-assembly requested but 'samples' is not an inline dict with 'rhizo' entries."
        )
    return [hifi_path_for_sample(s) for s in RHIZO_SAMPLES]

##################################################
# CLUSTER / SLURM RESOURCES (compatible, safer version)
##################################################

from functools import partial

def _build_slurm_partitions(cfg):
    """
    Build SLURM_PARTITIONS from config.
    Returns: [[name, min_mem(MB), max_mem(MB), min_threads, max_threads], ...]
    """
    parts = []
    try:
        for _, specs in cfg.get("slurm_partitions", {}).items():
            name = specs.get("name", "")
            # Memory in GB → MB (use 1000 to match your original behaviour)
            min_mem_mb = 1000 * int(specs.get("min_mem", 0))
            max_mem_mb = 1000 * int(specs.get("max_mem", 0))
            min_threads = int(specs.get("min_threads", 0))
            max_threads = int(specs.get("max_threads", 0))
            parts.append([name, min_mem_mb, max_mem_mb, min_threads, max_threads])
    except Exception:
        parts = []
    return parts

SLURM_PARTITIONS = _build_slurm_partitions(config)

def get_resource_real(
    wildcards,
    input,
    threads,
    attempt,
    SLURM_PARTITIONS="",
    mode="",
    mult=2,
    min_size=10000,
):
    """
    Drop-in replacement for original get_resource_real.
    Compatible with Snakemake dynamic resources and original behaviour.
    """
    def _return(mem, partition, thr, mode):
        if mode == "mem":
            return int(mem)
        if mode == "partition":
            return partition
        if mode == "threads":
            return int(thr)
        return int(mem)  # default

    # Estimate memory (MB)
    try:
        mem = max(
            (input.size // 1000000) * attempt * mult,
            attempt * min_size * mult
        )
    except Exception:
        mem = attempt * min_size * mult

    # No-cluster case
    if not SLURM_PARTITIONS or SLURM_PARTITIONS[0][0] == "":
        return _return(mem, "", threads, mode)

    # Filter partitions by capacity
    mem_ok = [p for p in SLURM_PARTITIONS if p[2] >= mem] or [max(SLURM_PARTITIONS, key=lambda x: x[2])]
    thr_ok = [p for p in mem_ok if p[4] >= threads] or [max(mem_ok, key=lambda x: x[4])]

    # Pick "tightest fit"
    name, min_mem, max_mem, min_thr, max_thr = min(
        thr_ok, key=lambda x: [x[1], x[3], x[2], x[4]]
    )

    # Clamp final values
    mem_final = min(max(mem, min_mem), max_mem)
    thr_final = min(max(threads, min_thr), max_thr)

    return _return(mem_final, name, thr_final, mode)

def get_resource(mode, **kwargs):
    return partial(
        get_resource_real,
        SLURM_PARTITIONS=SLURM_PARTITIONS,
        mode=mode,
        **kwargs
    )

resource_mem      = get_resource("mem")
resource_partition= get_resource("partition")
resource_threads  = get_resource("threads")

##################################################
# READ TYPE LOG (optional)
##################################################
# print("READ_TYPES:", READ_TYPES)
# if isinstance(cfg_samples, Mapping):
#     print("BULK_SAMPLES:", BULK_SAMPLES)
#     print("RHIZO_SAMPLES:", RHIZO_SAMPLES)
