"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-08-30]
Run: snakemake -s workflow/rules/metamdbg.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run metamdbg assembler on long reads
"""


############################################
rule lr_assembly:
    input:
        expand(os.path.join(RESULTS_DIR, "metamdbg/{sid}/{sid}.fasta"), sid=SAMPLES.index)
    output:
        touch("status/lr_assembly.done")


############################################
# localrules:


############################################
# Assembling the LR reads
rule metamdbg:
    input:
        lr=[lambda wildcards: SAMPLES.loc[wildcards.sid, "long_reads"]
    output:
        os.path.join(RESULTS_DIR, "metamdgb/{sid}/{sid}.fasta")
    conda:
        os.path.join(ENV_DIR, "metamdbg.yaml")
    threads:
        config['metamdbg']['threads']
    resources:
        mem_mb = resource_mem,
        slurm_partition = resource_partition
    log:
        os.path.join(RESULTS_DIR, "logs/metamdbg.{sid}.log")
    message:
        "Running metamdgb on {wildcards.sid}"
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    params:
        k=config['metamdbg']['k']
    benchmark:
        os.path.join(RESULTS_DIR, "benchmarks/metamdbg.{sid}.txt")
    shell:
        "(date && metaMDBG asm $(dirname {output}) {input} -t {threads} && "
        "metaMDBG gfa $(dirname {output}) {params.k} --contigpath --readpath && "
        "date) &> >(tee {log})"

