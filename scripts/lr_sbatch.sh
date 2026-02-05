#!/bin/bash -l

##############################
# SLURM
# NOTE: used for this script only, NOT for the snakemake call below

#SBATCH -J LR_decode
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=40-00:00:00
#SBATCH -p ei-long
#SBATCH --qos=qos-batch

##############################
# SNAKEMAKE

# conda env name
SMK_ENV="snakemake" # USER INPUT REQUIRED
# # number of cores for snakemake
# SMK_CORES=72
# number of jobs for snakemake
SMK_JOBS=16
# snakemake file
SMK_SMK="workflow/Snakefile"
# config file
SMK_CONFIG="config/lr_config.yaml" # USER INPUT REQUIRED
# slurm config file
SMK_SLURM="config/slurm.yaml"
# slurm profile file
SMK_PROFILE="profiles/slurm"
# # slurm cluster call
# SMK_CLUSTER="sbatch -p {cluster.partition} -q {cluster.qos} {cluster.explicit} -N {cluster.nodes} -n {cluster.n} -c {threads} -t {cluster.time} --job-name={cluster.job-name}"


##############################
# LAUNCHER

# activate the env
conda activate ${SMK_ENV}

# run the pipeline (without profile)
# snakemake -s ${SMK_SMK} --cores ${SMK_CORES} --local-cores 1 \
# --configfile ${SMK_CONFIG} --use-conda --conda-prefix ${CONDA_PREFIX}/pipeline \
# --cluster-config ${SMK_SLURM} --cluster "${SMK_CLUSTER}" --jobs "${SMK_JOBS}" --rerun-incomplete -rp

# run the pipeline (with profile)
snakemake --profile "${SMK_PROFILE}" --jobs "${SMK_JOBS}" --rerun-incomplete -rpk --unlock
snakemake --profile "${SMK_PROFILE}" --jobs "${SMK_JOBS}" --rerun-incomplete -rpk 

