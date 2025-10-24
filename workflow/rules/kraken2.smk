"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-08-30]
Run: snakemake -s workflow/rules/kraken2.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run Kraken2+BRACKEN on reads
"""


############################################
rule taxonomy:
    input:
        expand(os.path.join(RESULTS_DIR, "taxonomy/kraken2/{sid}_kraken.report"), sid=SAMPLES.index), 
        expand(os.path.join(RESULTS_DIR, "taxonomy/bracken/{sid}.bracken"), sid=SAMPLES.index),
        expand(os.path.join(RESULTS_DIR, "taxonomy/mpa_report/{sid}_mpa.tsv"), sid=SAMPLES.index),
        os.path.join(RESULTS_DIR, "taxonomy/mpa_report/combined_output.tsv"),
        os.path.join(RESULTS_DIR, "taxonomy/bracken/combined_bracken.txt")
    output:
        touch("status/taxonomy.done")


############################################
# localrules: phyloseq_input_kraken2


############################################
# Taxonomic classification using KRAKEN2
rule kraken2:
    input:
        dummy="status/databases.done",
        r1=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R1.fq"),
        r2=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R2.fq")
    output:
        report=os.path.join(RESULTS_DIR, "taxonomy/kraken2/{sid}_kraken.report"),
        summary=os.path.join(RESULTS_DIR, "taxonomy/kraken2/{sid}_kraken.out")
    conda:
        os.path.join(ENV_DIR, "kraken2.yaml")
    threads:
        config['kraken2']['threads']
    resources:
        mem_mb = resource_mem,
        slurm_partition = resource_partition
    params:
        db=config['kraken2']['db'],
        confidence=config['kraken2']['confidence']
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    log:
        os.path.join(RESULTS_DIR, "logs/kraken2/kraken2.{sid}.log")
    message:
        "Running kraken2 on {wildcards.sid}"
    shell:
        "(date && kraken2 --threads {threads} --db {params.db} --confidence {params.confidence} --paired --output {output.summary} --report {output.report} {input} && date) &> >(tee {log})"

# Running KRAKEN2+BRACKEN as suggested by KRAKEN2 website
rule bracken:
    input:
        report=rules.kraken2.output.report,
        r1=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R1.fq"),
        r2=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R2.fq")
    output:
        bracken=os.path.join(RESULTS_DIR, "taxonomy/bracken/{sid}.bracken"),
        report=os.path.join(RESULTS_DIR, "taxonomy/bracken/{sid}_bracken.report")
    threads:
        config['kraken2']['threads']
    resources:
        mem_mb = resource_mem,
        slurm_partition = resource_partition
    conda:
        os.path.join(ENV_DIR, "bracken.yaml")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    params:
        db=config['kraken2']['db'],
        read=config['kraken2']['read'],
        level=config['kraken2']['level'],
        bracken=config['bracken']['bin']
    log:
        os.path.join(RESULTS_DIR, "logs/bracken/bracken.{sid}.log")
    message:
        "Running kraken & bracken for {wildcards.sid}"
    shell:
        "(date && {params.bracken} -d {params.db} -i {input.report} -o {output.bracken} -w {output.report} -r {params.read} -l {params.level} && date)  &> >(tee {log})"

rule remove_uncultured:
    input:
        bracken=os.path.join(RESULTS_DIR, "taxonomy/bracken/{sid}.bracken")
    output:
        edited=os.path.join(RESULTS_DIR, "taxonomy/bracken/{sid}_edited.bracken")
    log:
        os.path.join(RESULTS_DIR, "logs/bracken/edited_bracken_{sid}")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Removing 'uncultured' taxa from bracken output from {wildcards.sid} due to combining issues"
    shell:
        "(date && grep -v 'uncultured' {input.bracken} | grep -v 'endosymbionts' | grep -v 'Incertae Sedis' > {output.edited} && date) &> >(tee {log})"

rule combine_bracken:
    input:
        bracken=expand(os.path.join(RESULTS_DIR, "taxonomy/bracken/{sid}_edited.bracken"), sid=SAMPLES.index)
    output:
        out=os.path.join(RESULTS_DIR, "taxonomy/bracken/combined_bracken.txt")
    conda:
        os.path.join(ENV_DIR, "python2.yaml")
    params:
        combine=config['bracken']['combine']
    resources:
        mem_mb = resource_mem,
        slurm_partition = resource_partition
    log:
        os.path.join(RESULTS_DIR, "logs/bracken/bracken_combine.log")
    message:
        "Combining all the output from BRACKEN"
    shell:
        "(date && python {params.combine} --files {input.bracken} -o {output.out} && date)  &> >(tee {log})"

#########################
### MPA-style report ###
rule mpa_report:
    input:
        report=os.path.join(RESULTS_DIR, "taxonomy/bracken/{sid}_bracken.report")
    output:
        mpa=os.path.join(RESULTS_DIR, "taxonomy/mpa_report/{sid}_mpa.tsv")
    conda:
        os.path.join(ENV_DIR, "bracken_new.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/bracken/mpa_{sid}.log")
    resources:
        mem_mb = resource_mem,
        slurm_partition = resource_partition
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Creating mpa-style report for {wildcards.sid}"
    shell:
        "(date && kreport2mpa.py -r {input.report} -o {output.mpa} && date)  &> >(tee {log})"

rule combine_mpa:
    input:
        mpa=expand(os.path.join(RESULTS_DIR, "taxonomy/mpa_report/{sid}_mpa.tsv"), sid=SAMPLES.index)
    output:
        combined=os.path.join(RESULTS_DIR, "taxonomy/mpa_report/combined_output.tsv")
    conda:
        os.path.join(ENV_DIR, "krakentools.yaml")
    resources:
        mem_mb = resource_mem,
        slurm_partition = resource_partition
    params:
        combine=os.path.join(SRC_DIR, "combine_mpa_modified.py")
    log:
        os.path.join(RESULTS_DIR, "logs/bracken/mpa_combine.log")
    message:
        "Creating a combined mpa-style report"
    shell:
        "(date && {params.combine} -i {input.mpa} -d $(dirname {output.combined}) && date)  &> >(tee {log})"


