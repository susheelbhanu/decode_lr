# -----------------------------------------------------------------
# SNAKEFILE - metaMDBG assemblies
# Using Singularity Container
# -----------------------------------------------------------------

import os

# --- Configuration & Paths ---
# Provide the absolute path to where you saved the .sif file
SIF_IMAGE = "/ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/decode_lr/singularity_cache/metamdbg_1.4.sif"

# --- Final targets ---
rule lr_assembly_all:
    input:
        expand(os.path.join(RESULTS_DIR, "assemblies_v1.4/hifi_individual/{sample}/contigs.fasta.gz"),
               sample=ALL_HIFI_SAMPLES),
        os.path.join(RESULTS_DIR, "assemblies_v1.4/hifi_coassembly_rhizo/contigs.fasta.gz"),
        os.path.join(RESULTS_DIR, "assemblies_v1.4/hifi_coassembly_bulk/contigs.fasta.gz"),
        os.path.join(RESULTS_DIR, "assemblies_v1.4/ont_individual_b26/contigs.fasta.gz"),
        os.path.join(RESULTS_DIR, "assemblies_v1.4/hybrid_b26/contigs.fasta.gz"),
        os.path.join(RESULTS_DIR, "assemblies_v1.4/hybrid_coassembly_bulk/contigs.fasta.gz")
    output:
        touch("status/lr_assembly_new_v1.4.done")

# --- Assembly Rules ---

rule metamdbg_hifi_individual:
    input:
        reads = get_hifi_reads
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies_v1.4/hifi_individual/{sample}/contigs.fasta.gz")
    params:
        opts = config["metamdbg"].get("opts", "")
    threads:
        config["metamdbg"]["threads"]
    singularity:
        SIF_IMAGE
    log:
        os.path.join(LOG_DIR, "metamdbg_hifi_individual/{sample}.log")
    shell:
        """
        (date && metaMDBG asm --in-hifi {input.reads} --out-dir $(dirname {output.contigs}) \
        --threads {threads} {params.opts} && date) &> >(tee {log})
        """

rule metamdbg_hifi_coassembly_rhizo:
    input:
        reads = get_rhizo_hifi_reads
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies_v1.4/hifi_coassembly_rhizo/contigs.fasta.gz")
    params:
        opts = config["metamdbg"].get("opts", "")
    threads:
        config["metamdbg"]["threads"]
    singularity:
        SIF_IMAGE
    log:
        os.path.join(LOG_DIR, "metamdbg_hifi_coassembly_rhizo.log")
    shell:
        """
        (date && metaMDBG asm --in-hifi {input.reads} --out-dir $(dirname {output.contigs}) \
        --threads {threads} {params.opts} && date) &> >(tee {log})
        """

rule metamdbg_hifi_coassembly_bulk:
    input:
        reads = get_bulk_hifi_reads
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies_v1.4/hifi_coassembly_bulk/contigs.fasta.gz")
    params:
        opts = config["metamdbg"].get("opts", "")
    threads:
        config["metamdbg"]["threads"]
    singularity:
        SIF_IMAGE
    log:
        os.path.join(LOG_DIR, "metamdbg_hifi_coassembly_bulk.log")
    shell:
        """
        (date && metaMDBG asm --in-hifi {input.reads} --out-dir $(dirname {output.contigs}) \
        --threads {threads} {params.opts} && date) &> >(tee {log})
        """

rule metamdbg_ont_individual:
    input:
        reads = ONT_FILE
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies_v1.4/ont_individual_b26/contigs.fasta.gz")
    params:
        opts = config["metamdbg"].get("opts", "")
    threads:
        config["metamdbg"]["threads"]
    singularity:
        SIF_IMAGE
    log:
        os.path.join(LOG_DIR, "metamdbg_ont_individual_b26.log")
    shell:
        """
        (date && metaMDBG asm --in-ont {input.reads} --out-dir $(dirname {output.contigs}) \
        --threads {threads} {params.opts} && date) &> >(tee {log})
        """

rule metamdbg_hybrid_b26:
    input:
        hifi = lambda wc: hifi_path_for_sample(HYBRID_INDIVIDUAL),
        ont  = ONT_FILE
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies_v1.4/hybrid_b26/contigs.fasta.gz")
    params:
        opts = config["metamdbg"].get("opts", "")
    threads:
        config["metamdbg"]["threads"]
    singularity:
        SIF_IMAGE
    log:
        os.path.join(LOG_DIR, "metamdbg_hybrid_b26.log")
    shell:
        """
        (date && metaMDBG asm --in-ont {input.hifi} {input.ont} \
        --out-dir $(dirname {output.contigs}) --threads {threads} {params.opts} && date) &> >(tee {log})
        """

rule metamdbg_hybrid_coassembly_bulk:
    input:
        ont_read = ONT_FILE
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies_v1.4/hybrid_coassembly_bulk/contigs.fasta.gz")
    params:
        opts          = config["metamdbg"].get("opts", ""),
        hifi_list_str = get_bulk_hifi_reads
    threads:
        config["metamdbg"]["threads"]
    singularity:
        SIF_IMAGE
    log:
        os.path.join(LOG_DIR, "metamdbg_hybrid_coassembly_bulk.log")
    shell:
        """
        (date && metaMDBG asm --in-ont {params.hifi_list_str} {input.ont_read} \
        --out-dir $(dirname {output.contigs}) --threads {threads} {params.opts} && date) &> >(tee {log})
        """
