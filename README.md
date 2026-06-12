# 🧬 DECODE-LR: Long-Read Metagenomics Workflow

## Overview
**DECODE-LR** is a modular [Snakemake](https://snakemake.readthedocs.io/) workflow for **metagenomic assembly, binning, annotation, and taxonomic profiling**, developed in the Molecular Ecology group at **UKCEH** and optimised for the **EI-HPC** (SLURM) cluster.

It is built around **long-read** data (**PacBio HiFi** and **Oxford Nanopore (ONT)**, including HiFi+ONT hybrid assemblies) but also supports short-read inputs through a sample sheet. The workflow is split into several independent entry points so each stage can be run on its own:

| Entry point | Snakefile | What it does |
|-------------|-----------|--------------|
| **Main** | `workflow/Snakefile` | Long-read assembly, binning, and phyloFlash (step-driven) |
| **Annotation** | `workflow/Snakefile_annotation` | LongFlow-style functional + taxonomic annotation of assemblies |
| **Standalone binning** | `workflow/rules/Snakefile_binning_standalone.smk` | Multi-binner MAG recovery + QC on a given assembly |
| **Benchmark** | `workflow/Snakefile_benchmark` | prodigal-gv vs pyrodigal-gv vs fetchMGs ORF/marker-gene comparison |

The pipeline is fully SLURM-aware (one profile per stage) and runs tools through **Singularity** containers and/or **conda** environments.

---

## 🔬 What the workflow does

**Assembly** (`rules/metamdbg_updated.smk`)
- Assembles HiFi, ONT, and hybrid (HiFi+ONT) reads with [metaMDBG **v1.4**](https://github.com/GaetanBenoitDev/metaMDBG) via a Singularity image.
- Modes: individual HiFi samples, group co-assemblies (`bulk` / `rhizo`), individual ONT, and hybrid assemblies.

**Binning** (`rules/lr_binning.smk`, `rules/Snakefile_binning_standalone.smk`)
- Read mapping with [minimap2](https://github.com/lh3/minimap2) (long reads) / [BWA](https://github.com/lh3/bwa) (short reads) and `samtools`.
- MAG recovery with [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/), [SemiBin2](https://github.com/BigDataBiology/SemiBin), and [CONCOCT](https://github.com/BinPro/CONCOCT).
- MAG QC and taxonomy with [CheckM2](https://github.com/chklovski/CheckM2), [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk), and FastTree.

**Annotation** (`rules/annotation.smk`) — adapted from **LongFlow**
- ORF calling with prodigal-gv, split across parallel batches.
- Single-copy gene (SCG) detection, ORF/SCG clustering (mmseqs2), and per-component quality.
- Functional annotation: DIAMOND vs CARD / IntI (ARGs, integrons), KO annotation via kofamscan.
- Taxonomy: CAT, geNomad (mobile/viral elements), and 16S rRNA classification via barrnap + BLCA.

**Taxonomic profiling** (`rules/phyloflash.smk`)
- [phyloFlash](https://github.com/HRGV/phyloFlash) SSU rRNA profiling per sample, with merged long- and wide-format abundance tables.

**MAG coverage** (`rules/mag_collection_cov_updated.smk`)
- Cross-sample MAG/dMAG coverage via ORF databases and [sylph](https://github.com/bluenote-1577/sylph) sketching/profiling.

**Benchmarking** (`rules/benchmark.smk`)
- Compares prodigal-gv vs pyrodigal-gv (whole fasta vs chunked) and rpsblast/hmmsearch vs [fetchMGs](https://github.com/motu-tool/fetchMGs) on a contig subset, writing Snakemake `benchmark` (time/memory) records.

---

## 📁 Repository Structure

```
decode_lr/
├── config/
│   ├── config.yaml                       # main (metag) config — used by workflow/Snakefile
│   ├── lr_config.yaml                    # long-read + annotation config
│   ├── binning_config.yaml               # standalone binning config
│   ├── mag_cov.yaml                       # MAG coverage config
│   ├── skyline_samples.tsv / samples.tsv # sample sheets (Sample_ID, sR1, sR2, long_reads, group)
│   └── rhizo_samples.txt
├── envs/                                 # conda environment YAMLs (metamdbg, mapping, metabat2, …)
├── schemas/                              # config.schema.yaml + samples.schema.yaml (validated at runtime)
├── workflow/
│   ├── Snakefile                         # main entry point (step-driven)
│   ├── Snakefile_annotation              # annotation entry point
│   ├── Snakefile_benchmark               # ORF-caller benchmark entry point
│   └── rules/
│       ├── init.smk                      # config validation, sample loading, resource helpers
│       ├── metamdbg_updated.smk          # metaMDBG v1.4 assemblies (Singularity)
│       ├── metamdbg.smk                  # legacy metaMDBG v1.3 (conda) assemblies
│       ├── lr_binning.smk                # mapping + MetaBAT2 binning
│       ├── Snakefile_binning_standalone.smk  # multi-binner + CheckM2 + GTDB-Tk
│       ├── annotation.smk                # LongFlow annotation rules
│       ├── benchmark.smk                 # prodigal vs pyrodigal vs fetchMGs
│       ├── phyloflash.smk                # phyloFlash SSU profiling
│       ├── mag_collection_cov_updated.smk    # MAG coverage via sylph
│       └── mag_collection_cov.smk        # legacy MAG coverage
├── profiles/
│   ├── slurm/                            # main workflow SLURM profile
│   ├── slurm_annotation/                 # annotation profile
│   ├── slurm_binning/                    # standalone binning profile
│   └── slurm_mag_cov/                    # MAG coverage profile
├── scripts/                              # SLURM launchers + helpers
│   ├── sbatch.sh                         # launches workflow/Snakefile
│   ├── lr_sbatch.sh                      # launches Snakefile_annotation
│   ├── binning_sbatch.sh                 # launches standalone binning
│   ├── benchmark.sh                      # launches Snakefile_benchmark
│   ├── sbatch_mag_cov.sh                 # launches MAG coverage
│   ├── metamdbg.sh, run_phyloflash.sh, prep_fastqs_by_tube.sh, extract_metaMDBG_stats.sh
└── LICENSE
```

> Outputs (assemblies, bins, annotation, logs, benchmarks) are written under the `results_dir` / `log_dir` / `benchmark_dir` paths set in the config — **not** inside the repo.

---

## ⚙️ Configuration

There are several config files; each entry point reads the relevant one. The **`steps`** list in the config controls which stages of the main workflow run.

### `config/config.yaml` (main / metag)
```yaml
# steps: ["profiles", "assembly", "mags", "analysis", "phyloflash"]
steps: ["profiles", "phyloflash"]

tag: "skyline"
workdir:     "/ei/projects/.../analysis"        # USER_INPUT
dbsdir:      "/ei/.project-scratch/.../databases"
samples:     "config/skyline_samples.tsv"        # USER_INPUT
results_dir: "/ei/projects/.../metag/results"
log_dir:     "/ei/projects/.../metag/logs"
benchmark_dir: "/ei/projects/.../metag/benchmarks"
```

### `config/lr_config.yaml` (long-read + annotation)
Recognised step keys: `lr_assembly`, `lr_mags`, `annotation`, `analysis`.
```yaml
steps: ["lr_assembly", "annotation"]

tag: "skyline_lr"
results_dir: "/ei/projects/.../lr_results_20260202"

# Inputs
hifi_dir: "/ei/projects/.../HiFi"
ont_file: "/ei/projects/.../ONT/CEHSoil_ONT_filtered_trimmed.fastq.gz"
hifi_file_pattern: "{sample}/{sample}_trimmed.fastq.gz"

# Sample groups
samples:
  bulk:  ["b10_d3_0200_con", "b20_d3_1400_bio", "b26_t3_con", ...]
  rhizo: ["r15_d3_1400_con", "r30_d3_0200_bio", ...]
hybrid_individual_sample: "b26_t3_con"

# Tools
metamdbg: { threads: 64, mem_gb: 750, k_min: 15, k_max: 27, opts: "" }
minimap2: { hifi_preset: "-ax map-hifi", ont_preset: "-ax map-ont" }
metabat2: { min_contig_length: 2000, opts: "--unbinned" }

# Annotation: per-assembly contig paths + database locations
longflow_dir:   "/ei/.../LongFlow"
annotation_dir: "/ei/.../annotation"
assembly_paths:
  hifi_individual_b26_t3_con: "/ei/.../assemblies/hifi_individual/b26_t3_con/contigs.fasta.gz"
  hybrid_b26:                 "/ei/.../assemblies/hybrid_b26/contigs.fasta.gz"
  # ... one flat key per assembly
annotation:
  prodigal_splits: 20
  diamond:  { CARD: {...}, IntI: {...} }   # leave a tool's paths empty to skip it
  kofamscan: {...}
  cat_path: "..."   ; cat_db: "..."
  genomad_db: "..."
  blca: {...}
```

> 💡 All paths must be **absolute** on EI-HPC. Settings tagged `# USER_INPUT` must be set before running, and user-specific paths should **not** be committed.

---

## 🛢️ Singularity Setup (EI-HPC)

Compute nodes (e.g. `ei-medium`) have **no internet access**, so containers must be pre-pulled on a login node into a shared cache, then reused via `--singularity-prefix`.

```bash
# 1. Create a shared cache directory
export SINGULARITY_CACHEDIR=/hpc-home/<user>/singularity_cache
mkdir -p "$SINGULARITY_CACHEDIR"
cd "$SINGULARITY_CACHEDIR"

# 2. Pull the containers used by the workflow (login node only)
singularity pull metamdbg_1.4.sif        docker://...metamdbg:1.4        # assembly
singularity pull minimap2_2.28.sif       docker://...minimap2:2.28       # mapping
singularity pull samtools_1.19.2.sif     docker://...samtools:1.19.2
singularity pull metabat2_2.15.sif       docker://...metabat2:2.15       # binning
singularity pull semibin2.sif            docker://...semibin            #   "
singularity pull concoct_1.1.0.sif       docker://...concoct:1.1.0      #   "
singularity pull checkm2.sif             docker://...checkm2            # MAG QC
singularity pull gtdbtk_2.4.0.sif        docker://...gtdbtk:2.4.0       #   "
singularity pull prodigal-gv_2.11.0.sif  docker://quay.io/biocontainers/prodigal-gv:2.11.0--h577a1d6_5  # annotation
singularity pull pyrodigalgv_0.3.2.sif   docker://...pyrodigal-gv:0.3.2 # benchmark
singularity pull fetchmgs.sif            docker://...fetchmgs           # benchmark
singularity pull mmseqs2_latest.sif      docker://...mmseqs2
singularity pull diamond_2.0.9.sif       docker://...diamond:2.0.9
singularity pull blast_2.16.0.sif        docker://quay.io/biocontainers/blast:2.16.0--hc155240_2
singularity pull barrnap_0.9.sif         docker://...barrnap:0.9
singularity pull genomad_1.8.0.sif       docker://...genomad:1.8.0
singularity pull phyloflash_3.3b1.sif    docker://...phyloflash:3.3b1
```

The profiles already enable Singularity (`use-singularity: True`) and point at the cache via `singularity-prefix`. Most rules also bind the cluster filesystems with `--singularity-args "--bind /ei,/hpc-home"`.

---

## 🧪 Conda Environments

Conda env YAMLs live in `envs/` (e.g. `metamdbg.yaml`, `mapping.yaml`, `metabat2.yaml`, `snakemake.yaml`, plus annotation/classification tools). The profiles set `use-conda: True` with a shared `conda-prefix`, so Snakemake builds envs on first run.

```bash
# Snakemake driver environment (used to launch the pipeline)
conda env create -f envs/snakemake.yaml
conda activate snakemake
```

---

## 🧰 SLURM Profiles

Each entry point has its own profile under `profiles/` (e.g. `profiles/slurm/config.yaml`). Resources are computed per-rule from the `slurm_partitions` block in the config, with profile-level fallbacks:

```yaml
default-resources:
  - slurm_account=tgac
  - slurm_partition="ei-long"
  - mem_mb=104857
  - time="1-00:00:00"
cluster: "sbatch --partition {resources.slurm_partition} --cpus-per-task={threads} --mem {resources.mem_mb} --qos {cluster.qos} --time {cluster.time} -J {cluster.job-name} -o slurm/%x-%j.out -e slurm/%x-%j.err"
use-conda: True
use-singularity: True
singularity-prefix: "/hpc-home/<user>/singularity_cache"   # USER_INPUT
```

> Set the `# USER_INPUT` fields (account, conda/singularity cache, `cluster-config` absolute path) before submitting.

---

## 🚀 Running the workflow

All stages are submitted with the launcher scripts in `scripts/`, which `conda activate` the Snakemake env and call Snakemake with the right Snakefile + profile.

```bash
# Main workflow (assembly / mags / phyloflash, per the `steps` list in config)
sbatch scripts/sbatch.sh

# Annotation of pre-computed assemblies
sbatch scripts/lr_sbatch.sh

# Standalone multi-binner MAG recovery + QC
sbatch scripts/binning_sbatch.sh

# MAG coverage profiling (sylph)
sbatch scripts/sbatch_mag_cov.sh

# ORF-caller / marker-gene benchmark
sbatch scripts/benchmark.sh
```

To run a stage manually (e.g. annotation), the underlying call looks like:

```bash
snakemake -s workflow/Snakefile_annotation \
    --profile profiles/slurm_annotation \
    --singularity-prefix /hpc-home/<user>/singularity_cache \
    --jobs 16 --use-singularity \
    --singularity-args "--bind /ei,/hpc-home" \
    --rerun-incomplete -rpk
```

Always **dry-run first** with `-n` to inspect the job graph and resource selection.

---

## 📂 Output Summary

| Step | Output (under `results_dir`) | Description |
|------|------------------------------|-------------|
| Assemblies | `assemblies_v1.4/hifi_individual/{sample}/contigs.fasta.gz` | metaMDBG v1.4 per-sample HiFi contigs |
| Co-assemblies | `assemblies_v1.4/hifi_coassembly_{bulk,rhizo}/` | Group co-assemblies |
| ONT / hybrid | `assemblies_v1.4/ont_individual_*`, `assemblies_v1.4/hybrid_*` | ONT and HiFi+ONT hybrid assemblies |
| Mapping | `mapping/{assembly_type}/` | Sorted, indexed BAMs |
| Binning | `binning/{run}/` | MetaBAT2 / SemiBin2 / CONCOCT bins + CheckM2 + GTDB-Tk |
| Annotation | `annotation/` | ORFs, SCGs, DIAMOND/CARD, KO, CAT, geNomad, 16S taxonomy |
| phyloFlash | `phyloflash/` + merged long/wide abundance tables | SSU rRNA community profiles |
| MAG coverage | sylph profiles / coverage matrices | Cross-sample MAG abundance |
| Benchmarks | `benchmark_dir/` | Snakemake time/memory records |

---

## 🧠 Tips
- Submit through the `scripts/*.sh` launchers — do **not** run Snakemake directly on the login node.
- Use `--rerun-incomplete` (already set in the profiles) to safely resume interrupted runs.
- Dry-run and inspect resources before launching:
  ```bash
  snakemake -s workflow/Snakefile --profile profiles/slurm -npr
  ```
- To skip an optional annotation tool, leave its database paths empty in `lr_config.yaml`.

---

## 📘 Citation

Please cite the core tools used in your analysis:

- **metaMDBG** — long-read metagenome assembly
- **minimap2** / **BWA** / **samtools** — read mapping
- **MetaBAT2**, **SemiBin2**, **CONCOCT** — binning
- **CheckM2**, **GTDB-Tk** — MAG QC and taxonomy
- **prodigal-gv** / **pyrodigal-gv**, **mmseqs2**, **DIAMOND**, **kofamscan**, **CAT**, **geNomad**, **fetchMGs**, **BLCA** — annotation
- **phyloFlash**, **sylph** — taxonomic profiling and coverage
- **Snakemake** — workflow management

---

*Author: Susheel Bhanu Busi — Molecular Ecology group, UKCEH.*
