# Initialization of a snakemake workflow
# Do not include here variables/settings which should/cannot be shared by all workflows

##################################################
# MODULES

import os
import re
import os.path
import pandas
from snakemake.utils import validate

##################################################
# CONFIG

# Config validation
validate(config, srcdir("../../schemas/config.schema.yaml"))
# Sample table (tab-separated, w/ header, 1st column is sample ID)
SAMPLES = pandas.read_csv(config["samples"], header=0, sep="\t").set_index("Sample_ID", drop=False)
# Sample table validation
validate(SAMPLES, srcdir("../../schemas/samples.schema.yaml"))

# Extract sample list for wildcards
SAMPLE = SAMPLES.index.tolist()

##################################################
# PATHS

SRC_DIR = srcdir("../../scripts") # add. scripts
ENV_DIR = srcdir("../../envs") # conda env. yaml files
MOD_DIR = srcdir("../../submodules") # git submodules
DBS_DIR = config["dbsdir"] # database folder

MDATA_DIR = os.path.abspath(config["metadata"]) # path to (meta)data

SAMPLES_FILE = os.path.abspath(config["samples"]) # samples: IDs and raw FASTQ files

OLDPWD = os.path.abspath(os.getcwd()) # PWD before changing the working directory

RESULTS_DIR = os.path.abspath(config["results_dir"])


# NEW PATHS
LOG_DIR = os.path.abspath(config["log_dir"])
BENCHMARK_DIR = os.path.abspath(config["benchmark_dir"])

BULK_SAMPLES = os.path.abspath(config["samples"]["bulk"])
RHIZO_SAMPLES = os.path.abspath(config["samples"]["rhizo"])
ALL_HIFI_SAMPLES = BULK_SAMPLES + RHIZO_SAMPLES

HYBRID_INDIVIDUAL = os.path.abspath(config["hybrid_individual_sample"])
ONT_FILE = os.path.abspath(config["ont_file"])

##################################################
# EXECUTION

# default executable for snakmake
shell.executable("bash")

# working directory
workdir:
    config["workdir"]

###################################################
# PARAMS

# File extensions of index files created by BWA
BWA_IDX_EXT = ["amb", "ann", "bwt", "pac", "sa"]

## Read types
#READ_TYPES = ["short"] # USER_INPUT: "short", "long"
## When `long reads` are specified`
#if SAMPLES["sR1"].notna() and SAMPLES["long_reads"].notna():
#    READ_TYPES = ["short", "long"]

READ_TYPES = []

# Check which columns exist
has_sR1_col = "sR1" in SAMPLES.columns
has_sR2_col = "sR2" in SAMPLES.columns
has_long_col = "long_reads" in SAMPLES.columns

# Detect short-read data (both R1 and R2 present for at least one sample)
has_short = (
    has_sR1_col and has_sR2_col
    and (SAMPLES["sR1"].notna() & SAMPLES["sR2"].notna()).any()
)

# Detect long-read data (any non-null entries)
has_long = has_long_col and SAMPLES["long_reads"].notna().any()

# Set read types based on what's present
if has_short:
    READ_TYPES.append("short")
if has_long:
    READ_TYPES.append("long")

# Safety check
if not READ_TYPES:
    raise WorkflowError(
        "No valid reads detected - check your sample sheet for sR1/sR2/long_reads columns."
    )


###################################################
# CLUSTER / SLURM RESOURCES (compatible, safer version)
# Drop-in replacement for the original get_resource_real

from functools import partial

# ----------------------------------------------------
# Build SLURM_PARTITIONS from config
# Format: [[name, min_mem(MB), max_mem(MB), min_threads, max_threads]]
# ----------------------------------------------------
def _build_slurm_partitions(cfg):
    parts = []
    try:
        for _, specs in cfg.get("slurm_partitions", {}).items():
            name = specs.get("name", "")
            # Original code: memory in GB, multiply by 1000 (not 1024)
            min_mem_mb = 1000 * int(specs.get("min_mem", 0))
            max_mem_mb = 1000 * int(specs.get("max_mem", 0))
            min_threads = int(specs.get("min_threads", 0))
            max_threads = int(specs.get("max_threads", 0))
            parts.append([name, min_mem_mb, max_mem_mb, min_threads, max_threads])
    except Exception:
        parts = []
    return parts

SLURM_PARTITIONS = _build_slurm_partitions(config)

# ----------------------------------------------------
# Main function (compatible with original get_resource_real)
# ----------------------------------------------------
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
    Compatible with Snakemake dynamic resources and old behaviour.
    """

    # Return helper
    def _return(mem, partition, thr, mode):
        if mode == "mem":
            return int(mem)
        if mode == "partition":
            return partition
        if mode == "threads":
            return int(thr)
        return int(mem)  # default fallback

    # Estimate memory (MB)
    try:
        mem = max(
            (input.size // 1000000) * attempt * mult,
            attempt * min_size * mult
        )
    except Exception:
        # Fall back if .size not available
        mem = attempt * min_size * mult

    # Handle no-cluster case (matches original)
    if not SLURM_PARTITIONS or SLURM_PARTITIONS[0][0] == "":
        return _return(mem, "", threads, mode)

    # --- Partition selection logic (same as original) ---
    # Step 1: filter partitions with enough max_mem
    mem_ok = [p for p in SLURM_PARTITIONS if p[2] >= mem]
    if not mem_ok:
        mem_ok = [max(SLURM_PARTITIONS, key=lambda x: x[2])]

    # Step 2: filter on threads
    thr_ok = [p for p in mem_ok if p[4] >= threads]
    if not thr_ok:
        thr_ok = [max(mem_ok, key=lambda x: x[4])]

    # Step 3: pick the "tightest fit"
    name, min_mem, max_mem, min_thr, max_thr = min(
        thr_ok, key=lambda x: [x[1], x[3], x[2], x[4]]
    )

    # Step 4: clamp final mem/threads
    mem_final = min(max(mem, min_mem), max_mem)
    thr_final = min(max(threads, min_thr), max_thr)

    return _return(mem_final, name, thr_final, mode)

# ----------------------------------------------------
# Partial wrappers (same usage as before)
# ----------------------------------------------------
def get_resource(mode, **kwargs):
    return partial(
        get_resource_real,
        SLURM_PARTITIONS=SLURM_PARTITIONS,
        mode=mode,
        **kwargs
    )

resource_mem = get_resource("mem")
resource_partition = get_resource("partition")
resource_threads = get_resource("threads")

# ----------------------------------------------------
# Functions for getting LR samples
# ----------------------------------------------------
def hifi_path_for_sample(sample: str) -> str:
    """
    Return the HiFi fastq path for a sample.
    Normal layout: {hifi_dir}/{sample}/{sample}_trimmed.fastq.gz
    Special-case for b26_t3_con: {hifi_dir}/{sample}/CEH_hifi_sample_trimmed.fastq.gz
    """
    base = os.path.abspath(config["hifi_dir"])
    p1 = os.path.join(base, sample, f"{sample}_trimmed.fastq.gz")
    p2 = os.path.join(base, sample, "CEH_hifi_sample_trimmed.fastq.gz")
    if os.path.exists(p1):
        return p1
    if os.path.exists(p2):
        return p2
    # Fall back to canonical path so Snakemake raises a clear MissingInputException
    return p1

def get_hifi_reads(wildcards):
    """Single-sample helper used by individual assembly rules."""
    return hifi_path_for_sample(wildcards.sample)

def get_bulk_hifi_reads(_wildcards=None):
    """Space-separated HiFi files for bulk co-assembly."""
    return " ".join(hifi_path_for_sample(s) for s in BULK_SAMPLES)

def get_rhizo_hifi_reads(_wildcards=None):
    """Space-separated HiFi files for rhizo co-assembly."""
    return " ".join(hifi_path_for_sample(s) for s in RHIZO_SAMPLES)

