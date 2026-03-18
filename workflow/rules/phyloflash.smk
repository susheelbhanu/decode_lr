from os.path import basename, dirname
import glob

# -----------------------------------------------------------------------------
# phyloFlash standalone rules
# -----------------------------------------------------------------------------

DATA = "/ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/metag/results"

PHYLOFLASH_DBHOME = "/ei/.project-scratch/0/0e51ef86-0156-4e79-ad12-c5411c0a5496/databases/phyloflash"
PHYLOFLASH_SIF = "/ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/decode_lr/singularity_cache/phyloflash_3.3b1.sif"

R1 = {
    basename(dirname(file)): file
    for file in glob.glob(f"{DATA}/preprocessed/reads/*/*_filtered.R1.fq")
}
R1.update({
    basename(dirname(file)): file
    for file in glob.glob(f"{DATA}/preprocessed/reads/*/*_filtered.R1.fq.gz")
})

R2 = {
    basename(dirname(file)): file
    for file in glob.glob(f"{DATA}/preprocessed/reads/*/*_filtered.R2.fq")
}
R2.update({
    basename(dirname(file)): file
    for file in glob.glob(f"{DATA}/preprocessed/reads/*/*_filtered.R2.fq.gz")
})

PHYLOFLASH_SAMPLES = sorted(set(R1) & set(R2))
PHYLOFLASH_DONE = f"{RESULTS_DIR}/phyloflash/{{sample}}/.phyloflash.done"

rule run_phyloflash:
    input:
        r1=lambda wc: R1[wc.sample],
        r2=lambda wc: R2[wc.sample]
    output:
        done=PHYLOFLASH_DONE
    threads: 24
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")
    params:
        dbhome=PHYLOFLASH_DBHOME,
        sif=PHYLOFLASH_SIF,
        lib=lambda wc: wc.sample,
        outdir=lambda wc: f"{RESULTS_DIR}/phyloflash/{wc.sample}"
    shell:
        r"""
        mkdir -p "{params.outdir}"

        singularity exec "{params.sif}" phyloFlash.pl \
          -lib "{params.lib}" \
          -read1 "{input.r1}" \
          -read2 "{input.r2}" \
          -CPUs {threads} \
          -dbhome "{params.dbhome}" \
          -almosteverything \
          -outdir "{params.outdir}"

        touch "{output.done}"
        """


PHYLOFLASH_TARGETS = expand(PHYLOFLASH_DONE, sample=PHYLOFLASH_SAMPLES)

rule all_phyloflash:
    input:
        PHYLOFLASH_TARGETS
