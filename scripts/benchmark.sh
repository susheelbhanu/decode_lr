#!/bin/bash -l

##############################
# SLURM  (orchestrator only — child jobs are submitted by the profile)
# NOTE: used for this script only, NOT for the snakemake child jobs below.
#
#SBATCH -J benchmark_ORFs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=2-00:00:00
#SBATCH -p ei-long
#SBATCH --qos=normal          # NOT qos-batch (doesn't exist here)
##############################

# ----------------------------------------------------------------------------
# PREREQUISITES (do these ONCE before sbatch):
#
# A) Fix the profile QOS — child jobs inherit qos from the cluster-config, and
#    profiles/slurm_annotation/slurm.yaml currently has __default__ qos: qos-batch
#    (invalid -> every job fails). Change that one line to:
#        qos: "normal"
#    (Partition is handled per-rule: the benchmark rules request ei-long, so the
#     profile's ei-largemem default is overridden automatically.)
#
# B) Pull the two NEW containers (login node, needs internet):
#    # pyrodigal-gv 0.3.2 (pure-python -> pyhdfd78af build):
#    curl -sL 'https://quay.io/api/v1/repository/biocontainers/pyrodigal-gv/tag/?limit=100&onlyActiveTags=true' \
#      | tr '{},' '\n' | grep '"name"' | grep '0\.3\.2'
#    cd /hpc-home/kar23heg/singularity_cache
#    singularity pull pyrodigalgv_0.3.2.sif \
#        https://depot.galaxyproject.org/singularity/pyrodigal-gv:0.3.2--pyhdfd78af_0
#
#    # fetchMGs:
#    curl -sL 'https://quay.io/api/v1/repository/biocontainers/fetchmgs/tag/?limit=100&onlyActiveTags=true' \
#      | tr '{},' '\n' | grep '"name"'
#    singularity pull fetchmgs.sif \
#        https://depot.galaxyproject.org/singularity/fetchmgs:<tag-from-above>
#    #   singularity exec fetchmgs.sif which fetchMGs.pl   # confirm entrypoint
# ----------------------------------------------------------------------------

##############################
# SNAKEMAKE
SMK_ENV="snakemake"                                   # USER INPUT if different
SMK_JOBS=12                                           # parallel child jobs
SMK_SMK="workflow/Snakefile_benchmark"
SMK_PROFILE="profiles/slurm_annotation"               # existing profile
SMK_SIF_CACHE="/hpc-home/kar23heg/singularity_cache"
##############################

# LAUNCHER
conda activate ${SMK_ENV}

mkdir -p slurm                                        # profile writes -o slurm/%x-%j.out

# container check
for sif in pyrodigalgv_0.3.2.sif fetchmgs.sif; do
    ls "${SMK_SIF_CACHE}/${sif}" >/dev/null 2>&1 || {
        echo "ERROR: ${SMK_SIF_CACHE}/${sif} not found — see PREREQUISITES B."; exit 1; }
done

# dry-run preview (uncomment):
# snakemake -s "${SMK_SMK}" --profile "${SMK_PROFILE}" \
#     --singularity-prefix "${SMK_SIF_CACHE}" --jobs "${SMK_JOBS}" -n

snakemake -s "${SMK_SMK}" --profile "${SMK_PROFILE}" --singularity-prefix "${SMK_SIF_CACHE}" --jobs "${SMK_JOBS}" --rerun-incomplete -rpk --unlock
snakemake -s "${SMK_SMK}" --profile "${SMK_PROFILE}" --singularity-prefix "${SMK_SIF_CACHE}" --jobs "${SMK_JOBS}" --rerun-incomplete -rpk 
