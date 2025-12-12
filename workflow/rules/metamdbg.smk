# -----------------------------------------------------------------
# SNAKEFILE - metaMDBG assemblies
# (handles special-case filename for b26_t3_con)
# -----------------------------------------------------------------

import os
import os.path

# --- Final targets ---
rule lr_assembly_all:
    input:
        # Individual HiFi assemblies
        expand(os.path.join(RESULTS_DIR, "assemblies/hifi_individual/{sample}/contigs.fasta.gz"),
               sample=ALL_HIFI_SAMPLES),
        # Co-assemblies
        os.path.join(RESULTS_DIR, "assemblies/hifi_coassembly_rhizo/contigs.fasta.gz"),
        os.path.join(RESULTS_DIR, "assemblies/hifi_coassembly_bulk/contigs.fasta.gz"),
        # ONT-only
        os.path.join(RESULTS_DIR, "assemblies/ont_individual_b26/contigs.fasta.gz"),
        # Hybrid
        os.path.join(RESULTS_DIR, "assemblies/hybrid_b26/contigs.fasta.gz"),
        os.path.join(RESULTS_DIR, "assemblies/hybrid_coassembly_bulk/contigs.fasta.gz")
    output:
        touch("status/lr_assembly.done")

# --- Assemblies (metaMDBG) ---
rule metamdbg_hifi_individual:
    """Assemble each HiFi sample individually."""
    input:
        reads = get_hifi_reads
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies/hifi_individual/{sample}/contigs.fasta.gz"),
        outdir  = directory(os.path.join(RESULTS_DIR, "assemblies/hifi_individual/{sample}"))
    params:
        mem  = config["metamdbg"]["mem_gb"],
        opts = config["metamdbg"].get("opts", "")
    threads:
        config["metamdbg"]["threads"]
    conda:
        "metamdbg" # os.path.join(ENV_DIR, "metamdbg.yaml")
    log:
        os.path.join(LOG_DIR, "metamdbg_hifi_individual/{sample}.log")
    benchmark:
        os.path.join(BENCHMARK_DIR, "metamdbg_hifi_individual/{sample}.txt")
    message:
        "metaMDBG asm: individual HiFi assembly of {wildcards.sample}"
    resources:
        mem_mb=resource_mem,
        slurm_partition=resource_partition
    shell:
        "(date && metaMDBG asm --in-hifi {input.reads} --out-dir {output.outdir} "
        "--threads {threads} {params.opts} && date) &> >(tee {log})"

rule metamdbg_hifi_coassembly_rhizo:
    """Co-assemble the rhizosphere HiFi samples."""
    input:
        reads = get_rhizo_hifi_reads
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies/hifi_coassembly_rhizo/contigs.fasta.gz"),
        outdir  = directory(os.path.join(RESULTS_DIR, "assemblies/hifi_coassembly_rhizo"))
    params:
        mem            = config["metamdbg"]["mem_gb"],
        opts           = config["metamdbg"].get("opts", ""),
        read_list_str  = get_rhizo_hifi_reads
    threads:
        config["metamdbg"]["threads"]
    conda:
        "metamdbg" # os.path.join(ENV_DIR, "metamdbg.yaml")
    log:
        os.path.join(LOG_DIR, "metamdbg_hifi_coassembly_rhizo.log")
    benchmark:
        os.path.join(BENCHMARK_DIR, "metamdbg_hifi_coassembly_rhizo.txt")
    message:
        "metaMDBG asm: HiFi co-assembly (rhizo)"
    resources:
        mem_mb=resource_mem,
        slurm_partition=resource_partition
    shell:
        "(date && metaMDBG asm --in-hifi {input.reads} --out-dir {output.outdir} "
        "--threads {threads} {params.opts} && date) &> >(tee {log})"

rule metamdbg_hifi_coassembly_bulk:
    """Co-assemble the bulk HiFi samples."""
    input:
        reads = get_bulk_hifi_reads
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies/hifi_coassembly_bulk/contigs.fasta.gz"),
        outdir  = directory(os.path.join(RESULTS_DIR, "assemblies/hifi_coassembly_bulk"))
    params:
        mem            = config["metamdbg"]["mem_gb"],
        opts           = config["metamdbg"].get("opts", ""),
        read_list_str  = get_bulk_hifi_reads
    threads:
        config["metamdbg"]["threads"]
    conda:
        "metamdbg" # os.path.join(ENV_DIR, "metamdbg.yaml")
    log:
        os.path.join(LOG_DIR, "metamdbg_hifi_coassembly_bulk.log")
    benchmark:
        os.path.join(BENCHMARK_DIR, "metamdbg_hifi_coassembly_bulk.txt")
    message:
        "metaMDBG asm: HiFi co-assembly (bulk)"
    resources:
        mem_mb=resource_mem,
        slurm_partition=resource_partition
    shell:
        "(date && metaMDBG asm --in-hifi {input.reads} --out-dir {output.outdir} "
        "--threads {threads} {params.opts} && date) &> >(tee {log})"

rule metamdbg_ont_individual:
    """Assemble the single ONT sample (b26) by itself."""
    input:
        reads = ONT_FILE
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies/ont_individual_b26/contigs.fasta.gz"),
        outdir  = directory(os.path.join(RESULTS_DIR, "assemblies/ont_individual_b26"))
    params:
        mem  = config["metamdbg"]["mem_gb"],
        opts = config["metamdbg"].get("opts", "")
    threads:
        config["metamdbg"]["threads"]
    conda:
        "metamdbg" # os.path.join(ENV_DIR, "metamdbg.yaml")
    log:
        os.path.join(LOG_DIR, "metamdbg_ont_individual_b26.log")
    benchmark:
        os.path.join(BENCHMARK_DIR, "metamdbg_ont_individual_b26.txt")
    message:
        "metaMDBG asm: ONT-only assembly (b26_t3_con)"
    resources:
        mem_mb=resource_mem,
        slurm_partition=resource_partition
    shell:
        "(date && metaMDBG asm --in-ont {input.reads} --out-dir {output.outdir} "
        "--threads {threads} {params.opts} && date) &> >(tee {log})"

rule metamdbg_hybrid_b26:
    """Hybrid assembly of b26_t3_con (HiFi + ONT)."""
    input:
        hifi = lambda wc: hifi_path_for_sample(HYBRID_INDIVIDUAL),
        ont  = ONT_FILE
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies/hybrid_b26/contigs.fasta.gz"),
        outdir  = directory(os.path.join(RESULTS_DIR, "assemblies/hybrid_b26"))
    params:
        mem  = config["metamdbg"]["mem_gb"],
        opts = config["metamdbg"].get("opts", "")
    threads:
        config["metamdbg"]["threads"]
    conda:
        "metamdbg" # os.path.join(ENV_DIR, "metamdbg.yaml")
    log:
        os.path.join(LOG_DIR, "metamdbg_hybrid_b26.log")
    benchmark:
        os.path.join(BENCHMARK_DIR, "metamdbg_hybrid_b26.txt")
    message:
        "metaMDBG asm: hybrid assembly (b26_t3_con: HiFi + ONT)"
    resources:
        mem_mb=resource_mem,
        slurm_partition=resource_partition
    shell:
        "(date && metaMDBG asm --in-ont {input.hifi} {input.ont} "
        "--out-dir {output.outdir} --threads {threads} {params.opts} && date) &> >(tee {log})"

rule metamdbg_hybrid_coassembly_bulk:
    """Hybrid co-assembly of bulk HiFi + ONT."""
    input:
        ont_read = ONT_FILE
    output:
        contigs = os.path.join(RESULTS_DIR, "assemblies/hybrid_coassembly_bulk/contigs.fasta.gz"),
        outdir  = directory(os.path.join(RESULTS_DIR, "assemblies/hybrid_coassembly_bulk"))
    params:
        mem          = config["metamdbg"]["mem_gb"],
        opts         = config["metamdbg"].get("opts", ""),
        hifi_list_str= get_bulk_hifi_reads
    threads:
        config["metamdbg"]["threads"]
    conda:
        "metamdbg" # os.path.join(ENV_DIR, "metamdbg.yaml")
    log:
        os.path.join(LOG_DIR, "metamdbg_hybrid_coassembly_bulk.log")
    benchmark:
        os.path.join(BENCHMARK_DIR, "metamdbg_hybrid_coassembly_bulk.txt")
    message:
        "metaMDBG asm: hybrid co-assembly (bulk: HiFi + ONT)"
    resources:
        mem_mb=resource_mem,
        slurm_partition=resource_partition
    shell:
        "(date && metaMDBG asm --in-ont {params.hifi_list_str} {input.ont_read} "
        "--out-dir {output.outdir} --threads {threads} {params.opts} && date) &> >(tee {log})"
