#!/bin/bash -l

##############################
# SLURM
#SBATCH -J phyloflash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=04-00:00:00
#SBATCH -p ei-medium
#SBATCH --qos=qos-batch

##############################
# SNAKEMAKE

SMK_ENV="snakemake"
SMK_JOBS=72
SMK_SMK="workflow/Snakefile"
SMK_CONFIG="config/config.yaml"
SMK_SLURM="config/slurm.yaml"
SMK_PROFILE="profiles/slurm"
SMK_TARGET="all_phyloflash"

##############################
# LAUNCHER

conda activate "${SMK_ENV}"

snakemake -s "${SMK_SMK}" \
  --configfile "${SMK_CONFIG}" \
  --profile "${SMK_PROFILE}" \
  --jobs "${SMK_JOBS}" \
  --rerun-incomplete -rpk \
  --unlock

snakemake -s "${SMK_SMK}" \
  --configfile "${SMK_CONFIG}" \
  --profile "${SMK_PROFILE}" \
  --jobs "${SMK_JOBS}" \
  --rerun-incomplete -rpk \
  "${SMK_TARGET}" all_phyloflash_html_abundance
