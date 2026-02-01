#!/bin/bash -l

##############################
# SLURM settings for launcher
##############################
#SBATCH -J BINNING_decode
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=20-00:00:00
#SBATCH -p ei-long
#SBATCH --qos=qos-batch


##############################
# Environment setup
##############################

conda activate snakemake
# module load apptainer   # or: module load singularity

##############################
# Directory binds for Singularity
##############################
# IMPORTANT: keep this as a single line so it passes cleanly via --singularity-args
BIND_DIRS="--bind $PWD --bind /ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7 --bind /ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7"
CACHE=/ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/decode_lr/singularity_cache
export SINGULARITY_CACHEDIR="$CACHE"
export APPTAINER_CACHEDIR="$CACHE"

##############################
# Run Snakemake
##############################

cd /ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/decode_lr || exit 1

SMK_SMK="workflow/rules/Snakefile_binning_standalone.smk"
SMK_PROFILE="profiles/slurm_binning"
SMK_CONFIG="config/binning_config.yaml"
SMK_JOBS=24


# Unlock (in case a previous run left a lock)
snakemake -s "${SMK_SMK}" --profile "${SMK_PROFILE}" --configfile "${SMK_CONFIG}" \
  --jobs "${SMK_JOBS}" --rerun-incomplete -rpk --unlock \
  --use-singularity --singularity-prefix "$CACHE" --singularity-args "${BIND_DIRS}"

# Run
snakemake -s "${SMK_SMK}" --profile "${SMK_PROFILE}" --configfile "${SMK_CONFIG}" \
  --jobs "${SMK_JOBS}" --rerun-incomplete -rpk \
  --use-singularity --singularity-prefix "$CACHE" --singularity-args "${BIND_DIRS}"
