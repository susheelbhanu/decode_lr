"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-08-30]
Run: snakemake -s workflow/rules/metaphlan.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run MetaPhlAn (v4.2.3) on paired-end reads and produce per-sample profiles + merged table
"""

############################################
# Final taxonomy target rule
rule metaphlan_all:
    input:
        expand(os.path.join(RESULTS_DIR, "taxonomy/metaphlan/{sid}_mpa.tsv"), sid=SAMPLES.index),
        os.path.join(RESULTS_DIR, "taxonomy/mpa_report/combined_metaphlan.tsv")
    output:
        touch("status/metaphlan.done")


############################################
localrules: metaphlan_merge


############################################
# Taxonomic classification using MetaPhlAn v4.2.3
rule metaphlan:
    input:
        r1 = os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R1.fq"),
        r2 = os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R2.fq")
    output:
        prof = os.path.join(RESULTS_DIR, "taxonomy/metaphlan/{sid}_mpa.tsv"),
        bt2  = os.path.join(RESULTS_DIR, "taxonomy/metaphlan/{sid}.bowtie2.bz2")
    conda:
        "metaphlan" #os.path.join(ENV_DIR, "metaphlan.yaml")
    threads:
        config["metaphlan"]["threads"]
    resources:
        mem_mb = resource_mem,
        slurm_partition = resource_partition
    params:
        bowtie2db = config["metaphlan"]["bowtie2db"],
        index     = config["metaphlan"]["index"],
        unknown   = "--unknown_estimation" if config["metaphlan"].get("unknown_estimation", False) else "",
        add_strain= "--add_merged_nreads --min_cu_len 2000 --stat_q 0.2" if config["metaphlan"].get("sensible_defaults", True) else ""
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    log:
        os.path.join(RESULTS_DIR, "logs/metaphlan/metaphlan.{sid}.log")
    message:
        "Running MetaPhlAn v4.2.3 on {wildcards.sid}"
    shell:
        "(date && metaphlan {input.r1},{input.r2} "
        "--input_type fastq "
        "--bowtie2db {params.bowtie2db} "
        "--index {params.index} "
        "--nproc {threads} "
        "--bowtie2out {output.bt2} "
        "{params.unknown} {params.add_strain} "
        "-o {output.prof} && date) &> >(tee {log})"


############################################
# Merge MetaPhlAn outputs into a single table
rule metaphlan_merge:
    input:
        expand(os.path.join(RESULTS_DIR, "taxonomy/metaphlan/{sid}_mpa.tsv"), sid=SAMPLES.index)
    output:
        merged = os.path.join(RESULTS_DIR, "taxonomy/mpa_report/combined_metaphlan.tsv")
    conda:
        "metaphlan" # os.path.join(ENV_DIR, "metaphlan.yaml")
    threads: 4
    resources:
        mem_mb = resource_mem,
        slurm_partition = resource_partition
    log:
        os.path.join(RESULTS_DIR, "logs/metaphlan/metaphlan_merge.log")
    message:
        "Merging MetaPhlAn profiles"
    shell:
        "(date && merge_metaphlan_tables.py {input} > {output.merged} && date) &> >(tee {log})"
