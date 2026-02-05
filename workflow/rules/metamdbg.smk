# -----------------------------------------------------------------
# SNAKEFILE - metaMDBG assemblies
# Using conda 'metamdbg v1.3' and manual executable path
# -----------------------------------------------------------------

import os

# --- Configuration & Paths ---
METAMDBG_EXE = "/hpc-home/kar23heg/tools/metaMDBG/build/bin/metaMDBG"
CONDA_ENV = "metamdbg1.3"

# --- Final targets ---
rule lr_assembly_all:
    input:
        expand(os.path.join(RESULTS_DIR, "assemblies/hifi_individual/{sample}/contigs.fasta.gz"),
               sample=ALL_HIFI_SAMPLES),
        os.path.join(RESULTS_DIR, "assemblies/hifi_coassembly_rhizo/contigs.fasta.gz"),
        os.path.join(RESULTS_DIR, "assemblies/hifi_coassembly_bulk/contigs.fasta.gz"),
        os.path.join(RESULTS_DIR, "assemblies/ont_individual_b26/contigs.fasta.gz"),
        os.path.join(RESULTS_DIR, "assemblies/hybrid_b26/contigs.fasta.gz"),
        os.path.join(RESULTS_DIR, "assemblies/hybrid_coassembly_bulk/contigs.fasta.gz")
    output:
        touch("status/lr_assembly_new.done")

# --- Assembly Rules ---

rule metamdbg_hifi_individual:
    input:
        reads = get_hifi_reads
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies/hifi_individual/{sample}/contigs.fasta.gz"),
        outdir  = directory(os.path.join(RESULTS_DIR, "assemblies/hifi_individual/{sample}"))
    params:
        opts = config["metamdbg"].get("opts", "")
    threads:
        config["metamdbg"]["threads"]
    conda:
        CONDA_ENV
    log:
        os.path.join(LOG_DIR, "metamdbg_hifi_individual/{sample}.log")
    shell:
        """
        # Ensure conda-installed time/minimap2 are prioritized
        export PATH=$CONDA_PREFIX/bin:$PATH

        (date && {METAMDBG_EXE} asm --in-hifi {input.reads} --out-dir {output.outdir} \
        --threads {threads} {params.opts} && date) &> >(tee {log})
        """

rule metamdbg_hifi_coassembly_rhizo:
    input:
        reads = get_rhizo_hifi_reads
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies/hifi_coassembly_rhizo/contigs.fasta.gz"),
        outdir  = directory(os.path.join(RESULTS_DIR, "assemblies/hifi_coassembly_rhizo"))
    params:
        opts = config["metamdbg"].get("opts", "")
    threads:
        config["metamdbg"]["threads"]
    conda:
        CONDA_ENV
    log:
        os.path.join(LOG_DIR, "metamdbg_hifi_coassembly_rhizo.log")
    shell:
        """
        export PATH=$CONDA_PREFIX/bin:$PATH
        (date && {METAMDBG_EXE} asm --in-hifi {input.reads} --out-dir {output.outdir} \
        --threads {threads} {params.opts} && date) &> >(tee {log})
        """

rule metamdbg_hifi_coassembly_bulk:
    input:
        reads = get_bulk_hifi_reads
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies/hifi_coassembly_bulk/contigs.fasta.gz"),
        outdir  = directory(os.path.join(RESULTS_DIR, "assemblies/hifi_coassembly_bulk"))
    params:
        opts = config["metamdbg"].get("opts", "")
    threads:
        config["metamdbg"]["threads"]
    conda:
        CONDA_ENV
    log:
        os.path.join(LOG_DIR, "metamdbg_hifi_coassembly_bulk.log")
    shell:
        """
        export PATH=$CONDA_PREFIX/bin:$PATH
        (date && {METAMDBG_EXE} asm --in-hifi {input.reads} --out-dir {output.outdir} \
        --threads {threads} {params.opts} && date) &> >(tee {log})
        """

rule metamdbg_ont_individual:
    input:
        reads = ONT_FILE
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies/ont_individual_b26/contigs.fasta.gz"),
        outdir  = directory(os.path.join(RESULTS_DIR, "assemblies/ont_individual_b26"))
    params:
        opts = config["metamdbg"].get("opts", "")
    threads:
        config["metamdbg"]["threads"]
    conda:
        CONDA_ENV
    log:
        os.path.join(LOG_DIR, "metamdbg_ont_individual_b26.log")
    shell:
        """
        export PATH=$CONDA_PREFIX/bin:$PATH
        (date && {METAMDBG_EXE} asm --in-ont {input.reads} --out-dir {output.outdir} \
        --threads {threads} {params.opts} && date) &> >(tee {log})
        """

rule metamdbg_hybrid_b26:
    input:
        hifi = lambda wc: hifi_path_for_sample(HYBRID_INDIVIDUAL),
        ont  = ONT_FILE
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies/hybrid_b26/contigs.fasta.gz"),
        outdir  = directory(os.path.join(RESULTS_DIR, "assemblies/hybrid_b26"))
    params:
        opts = config["metamdbg"].get("opts", "")
    threads:
        config["metamdbg"]["threads"]
    conda:
        CONDA_ENV
    log:
        os.path.join(LOG_DIR, "metamdbg_hybrid_b26.log")
    shell:
        """
        export PATH=$CONDA_PREFIX/bin:$PATH
        (date && {METAMDBG_EXE} asm --in-ont {input.hifi} {input.ont} \
        --out-dir {output.outdir} --threads {threads} {params.opts} && date) &> >(tee {log})
        """

rule metamdbg_hybrid_coassembly_bulk:
    input:
        ont_read = ONT_FILE
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies/hybrid_coassembly_bulk/contigs.fasta.gz"),
        outdir  = directory(os.path.join(RESULTS_DIR, "assemblies/hybrid_coassembly_bulk"))
    params:
        opts          = config["metamdbg"].get("opts", ""),
        hifi_list_str = get_bulk_hifi_reads
    threads:
        config["metamdbg"]["threads"]
    conda:
        CONDA_ENV
    log:
        os.path.join(LOG_DIR, "metamdbg_hybrid_coassembly_bulk.log")
    shell:
        """
        export PATH=$CONDA_PREFIX/bin:$PATH
        (date && {METAMDBG_EXE} asm --in-ont {params.hifi_list_str} {input.ont_read} \
        --out-dir {output.outdir} --threads {threads} {params.opts} && date) &> >(tee {log})
        """
