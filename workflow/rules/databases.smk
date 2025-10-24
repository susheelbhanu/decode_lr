"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-10-15]
Run: snakemake -s workflow/rules/databases.smk --configfile config/config.yaml --use-conda --cores 4 -rp
Latest modification:
Purpose: Download and unpack relevant databases
"""

############################################
rule databases:
    input:
        os.path.join(DBS_DIR, "kraken2/hash.k2d"),
        os.path.join(DBS_DIR, "kaiju/kaiju_db_nr_euk.fmi")
    output:
        touch("status/databases.done")


############################################
localrules: kraken_db


############################################
# KRAKEN database
rule kraken_db:
    input:
        url=config["kraken2"]["url"]
    output:
        os.path.join(DBS_DIR, "kraken2/hash.k2d")
    log:
        os.path.join(RESULTS_DIR, "logs/db_download/kraken.log")
#    resources:
#        mem_mb = resource_mem,
#        slurm_partition = resource_partition
    message:
        "KRAKEN2: database download"
    shell:
        "(date && "
        "cd $(dirname {output}) && wget {input.url} && "
        "tar -xzvf $(basename {input.url}) && "
        "date) &> {log}"

# KAIJU database
rule kaiju_db:
    input:
        url=config["kaiju"]["url"]
    output:
        os.path.join(DBS_DIR, "kaiju/kaiju_db_nr_euk.fmi")
    log:
        os.path.join(RESULTS_DIR, "logs/db_download/kaiju.log")
#    resources:
#        mem_mb = resource_mem,
#        slurm_partition = resource_partition
    message:
        "KAIJU: databse download"
    shell:
        "(date && "
        "cd $(dirname {output}) && wget {input.url} && "
        "tar -xzvf $(basename {input.url}) && "
        "date) &> {log}"


