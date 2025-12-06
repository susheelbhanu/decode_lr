# 🧬 DECODE-LR: Long-Read Assembly and Binning Workflow

## Overview
**DECODE-LR** is a Snakemake-based workflow for **long-read metagenomic assembly and binning**, optimised for **EI-HPC**.  
It supports **PacBio HiFi** and **Oxford Nanopore (ONT)** data, and automates:

1. **Assembly** using [metaMDBG](https://github.com/Roy-Lab/metaMDBG)  
2. **Read mapping** using [minimap2](https://github.com/lh3/minimap2)  
3. **Binning** using [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/)

The workflow is fully SLURM-aware and designed for scalable hybrid assemblies (HiFi + ONT).

---

## 📁 Repository Structure

```
decode_lr/
├── config/
│   └── lr_config.yaml                 # user-editable config
├── envs/                              # conda environment YAMLs
├── schemas/                           # JSON/YAML validation schemas
├── workflow/
│   ├── metamdbg_Snakefile             # metaMDBG assemblies
│   ├── metamdbg_binning_Snakefile     # mapping + binning
│   ├── rules/init.smk                 # setup, helpers, and resource logic
│   └── profiles/slurm/                # EI-HPC SLURM runtime configs
├── scripts/                           # helper scripts (if any)
├── lr_sbatch.sh                       # EI-HPC submission script
└── results/                           # all outputs (assemblies, bins, logs)
```

---

## ⚙️ Configuration (`config/lr_config.yaml`)

Example:

```yaml
# Input directories
hifi_dir: "/ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/HiFi"
ont_file: "/ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/ONT/CEHSoil_filtered_trimmed.fastq.gz"

# Sample groups
samples:
  bulk:
    - "b10_d3_0200_con"
    - "b20_d3_1400_bio"
    - "b26_t3_con"
  rhizo:
    - "r15_d3_1400_con"
    - "r30_d3_0200_bio"

# Hybrid assembly
hybrid_individual_sample: "b26_t3_con"

# Tool settings
metamdbg:
  threads: 48
  mem_gb: 250
  opts: ""

minimap2:
  threads: 16
  hifi_preset: "-ax map-hifi"
  ont_preset: "-ax map-ont"

metabat2:
  min_contig_length: 2000
  opts: "--unbinned"

# Optional SLURM partitions
slurm_partitions:
  highmem:
    name: "ei-largemem"
    min_mem: 750
    max_mem: 6000
```

💡 *Paths must be absolute on EI-HPC.*

---

# 🛢️ Singularity Setup (EI-HPC Software node)

Compute nodes such as `ei-medium` **do not have internet access**, so Singularity must use:

1. A **local project-owned cache**
2. **Local .sif images** pulled on a login node

---

## Create a Singularity cache directory

```bash
export SINGULARITY_CACHEDIR=/ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/decode_lr/singularity_cache
mkdir -p "$SINGULARITY_CACHEDIR"
```

## Create the containers folder

```bash
mkdir -p /ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/decode_lr/containers
cd /ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/decode_lr/containers
```

## Pull workflow containers (login node only)

```bash
singularity pull bwasamtools_1.10.sif      docker://quay.io/annacprice/bwasamtools:1.10
singularity pull bedtools_2.29.2.sif       docker://quay.io/annacprice/bedtools:2.29.2
singularity pull prodigal-gv_2.11.0.sif    docker://quay.io/biocontainers/prodigal-gv:2.11.0--h577a1d6_5
singularity pull blast_2.16.0.sif          docker://quay.io/biocontainers/blast:2.16.0--hc155240_2
singularity pull pythonenv_3.9.sif         docker://quay.io/annacprice/pythonenv:3.9
singularity pull sylph_0.9.0.sif           docker://quay.io/biocontainers/sylph:0.9.0--ha6fb395_0
```

## Use containers in Snakemake

```python
CONTAINER_DIR = "/ei/projects/.../decode_lr/containers"

BWASAMTOOLS_IMG = f"{CONTAINER_DIR}/bwasamtools_1.10.sif"
BEDTOOLS_IMG    = f"{CONTAINER_DIR}/bedtools_2.29.2.sif"
PRODIGAL_IMG    = f"{CONTAINER_DIR}/prodigal-gv_2.11.0.sif"
BLAST_IMG       = f"{CONTAINER_DIR}/blast_2.16.0.sif"
PYTHONENV_IMG   = f"{CONTAINER_DIR}/pythonenv_3.9.sif"
SYLPH_IMG       = f"{CONTAINER_DIR}/sylph_0.9.0.sif"
```

Then within rules:

```python
singularity: BWASAMTOOLS_IMG
```

---

## 🧪 Conda Environments

Create conda environments (once):

```bash
conda env create -f envs/metamdbg.yaml
conda env create -f envs/mapping.yaml
conda env create -f envs/metabat2.yaml
```

Activate the Snakemake environment before submission:

```bash
conda activate decode
```

---

## 🧰 EI-HPC Profile

`workflow/profiles/slurm/config.yaml` defines default resources and SLURM options, e.g.:

```yaml
default-resources:
  - slurm_partition="ei-medium,ei-long"
  - mem_mb=10000
  - time="6:00:00"
cluster: "sbatch --partition {resources.slurm_partition} --cpus-per-task={threads} --mem {resources.mem_mb} --qos {cluster.qos} --time {cluster.time} -J {cluster.job-name} -o slurm/%x-%j.out -e slurm/%x-%j.err"
```

Each rule computes:
```python
resources:
    mem_mb = resource_mem,
    slurm_partition = resource_partition
```

Fallbacks apply automatically from the profile if partitions are not explicitly set.

---

## 🚀 Execution on EI-HPC

The pipeline is launched through the provided SLURM wrapper script.

### 1️⃣ Submit the long-read assembly workflow
```bash
sbatch lr_sbatch.sh workflow/metamdbg_Snakefile config/lr_config.yaml
```

### 2️⃣ Submit the mapping + binning workflow
```bash
sbatch lr_sbatch.sh workflow/metamdbg_binning_Snakefile config/lr_config.yaml
```

### Example `lr_sbatch.sh`

```bash
#!/bin/bash
#SBATCH -p ei-medium
#SBATCH -t 72:00:00
#SBATCH -J decode_lr
#SBATCH -o slurm/%x-%j.out
#SBATCH -e slurm/%x-%j.err
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --qos=ei-medium

module load mambaforge/23.3.1
conda activate decode

snakemake -s $1 --configfile $2 --profile workflow/profiles/slurm --use-conda --cores 48
```

---

## 📂 Output Summary

| Step | Output directory | Description |
|------|------------------|--------------|
| Individual assemblies | `results/assemblies/hifi_individual/{sample}/` | metaMDBG contigs per HiFi sample |
| Co-assemblies (bulk/rhizo) | `results/assemblies/hifi_coassembly_*` | Combined assemblies by group |
| Hybrid assembly (HiFi + ONT) | `results/assemblies/hybrid_*` | Cross-technology assembly |
| Mapping | `results/mapping/{assembly_type}/` | BAMs, indices |
| Binning | `results/binning/{assembly_type}/` | MAG bins + completeness |

---

## 🧠 Tips
- Always submit via `sbatch lr_sbatch.sh` — do **not** run directly on the login node.  
- Use `--rerun-incomplete` to resume interrupted runs safely.  
- Inspect resource selection before launch:  
  ```bash
  snakemake -s workflow/metamdbg_Snakefile --configfile config/lr_config.yaml -npr --scheduler greedy
  ```

---

## 📘 Citation

Please cite the core tools used:

- **metaMDBG**
- **MetaBAT2**  
- **Minimap2**  
- **Snakemake**  

