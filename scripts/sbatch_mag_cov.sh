#!/bin/bash -l

##############################
# SLURM settings for launcher
##############################

#SBATCH -J MAG_cov
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --time=20-00:00:00
#SBATCH -p ei-largemem
#SBATCH --qos=qos-batch
#SBATCH --mem=600G 

##############################
# Environment setup
##############################

# Activate your Snakemake environment
conda activate snakemake

## Load Singularity / Apptainer module if required by your cluster
# module load apptainer   # or: module load singularity

##############################
# Directory binds for Singularity
##############################

BIND_DIRS="--bind $PWD \
           --bind /ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7 \
           --bind /ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7"

##############################
# Run Snakemake
##############################
cd /ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/decode_lr

snakemake -s workflow/rules/mag_collection_cov_updated.smk \
    --cores 92 --jobs 8 \
    --use-conda \
    --use-singularity \
    --singularity-args "$BIND_DIRS" \
    --rerun-incomplete -k \
    --unlock

snakemake -s workflow/rules/mag_collection_cov_updated.smk \
     --cores 92 --jobs 8 \
     --use-conda \
     --use-singularity \
     --singularity-args "$BIND_DIRS" \
     --rerun-incomplete -k 
