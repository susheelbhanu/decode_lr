#!/bin/bash -l

##############################
# SLURM
# NOTE: these resources are for the *launcher* only.
# Actual rule jobs are submitted via the Snakemake SLURM profile.
##############################

#SBATCH -J MAG_cov
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=20-00:00:00
#SBATCH -p ei-long
#SBATCH --mem=8G

##############################
# SNAKEMAKE
##############################

# conda env name
SMK_ENV="snakemake"

# number of concurrent cluster jobs Snakemake can submit
SMK_JOBS=8

# snakemake rule file
SMK_SMK="workflow/rules/mag_collection_cov_updated.smk"

# snakemake config file for this run
SMK_CONFIG="config/mag_cov.yaml"

# slurm profile
SMK_PROFILE="profiles/slurm_mag_cov"

##############################
# SINGULARITY / APPTAINER BINDS
##############################

BIND_DIRS="--bind $PWD \
           --bind /ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7 \
           --bind /ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7"

##############################
# LAUNCHER
##############################

# activate environment
conda activate "${SMK_ENV}"

# move to project directory
cd /ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/decode_lr || exit 1

# ensure local tmp exists (profile exports TMPDIR to jobs; rules/tools may use it)
mkdir -p /ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7/tmp

# unlock first (clears stale locks)
snakemake -s "${SMK_SMK}" \
    --profile "${SMK_PROFILE}" \
    --configfile "${SMK_CONFIG}" \
    --jobs "${SMK_JOBS}" \
    --use-conda \
    --use-singularity \
    --singularity-args "$BIND_DIRS" \
    --rerun-incomplete -rpk --unlock

# actual run
snakemake -s "${SMK_SMK}" \
    --profile "${SMK_PROFILE}" \
    --configfile "${SMK_CONFIG}" \
    --jobs "${SMK_JOBS}" \
    --use-conda \
    --use-singularity \
    --singularity-args "$BIND_DIRS" \
    --rerun-incomplete -rpk 
