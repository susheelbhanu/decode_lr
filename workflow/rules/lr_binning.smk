# -----------------------------------------------------------------
# SNAKEFILE - metaMDBG binning
# (handles special-case filename for b26_t3_con)
# -----------------------------------------------------------------

import os
import os.path

rule binning_all:
    input:
        # Hybrid + binning completion flags (as per your previous plan)
        os.path.join(RESULTS_DIR, "binning/hybrid_b26/bins.DONE"),
        os.path.join(RESULTS_DIR, "binning/hybrid_coassembly_bulk/bins.DONE")
    output:
        touch("status/lr_binning.done")

# --- Mapping (minimap2 + samtools) ---

def mapping_dir(assembly_type):
    return os.path.join(RESULTS_DIR, "mapping", assembly_type)

def binning_dir(assembly_type):
    return os.path.join(RESULTS_DIR, "binning", assembly_type)

rule minimap2_index:
    """Index assemblies for mapping."""
    input:
        os.path.join(RESULTS_DIR, "assemblies/{assembly_type}/contigs.fasta")
    output:
        os.path.join(mapping_dir("{assembly_type}"), "contigs.mmi")
    conda:
        os.path.join(ENV_DIR, "mapping.yaml")
    log:
        os.path.join(LOG_DIR, "minimap2_index/{assembly_type}.log")
    benchmark:
        os.path.join(BENCHMARK_DIR, "minimap2_index/{assembly_type}.txt")
    message:
        "minimap2: index {wildcards.assembly_type}"
    shell:
        "(date && minimap2 -d {output} {input} && date) &> >(tee {log})"

rule map_reads_to_assembly:
    input:
        idx      = os.path.join(mapping_dir("{assembly_type}"), "contigs.mmi"),
        assembly = os.path.join(RESULTS_DIR, "assemblies/{assembly_type}/contigs.fasta"),
        reads    = lambda wc: hifi_path_for_sample(wc.sample) if wc.read_tech == "hifi" else ONT_FILE
    output:
        # sample cannot contain '.', read_tech must be 'hifi' or 'ont'
        bam = os.path.join(mapping_dir("{assembly_type}"),
                           "{sample,[^\\.]+}.{read_tech,(hifi|ont)}.bam")
    params:
        preset = lambda wc: config["minimap2"]["hifi_preset"] if wc.read_tech == "hifi" else config["minimap2"]["ont_preset"]
    threads:
        config["minimap2"]["threads"]
    conda:
        os.path.join(ENV_DIR, "mapping.yaml")
    log:
        os.path.join(LOG_DIR, "map_reads/{assembly_type}/{sample}.{read_tech}.log")
    message:
        "minimap2: map {wildcards.read_tech} reads from {wildcards.sample} to {wildcards.assembly_type}"
    shell:
        "(date && minimap2 -t {threads} {params.preset} {input.assembly} {input.reads} | "
        "samtools view -@ {threads} -bS - > {output.bam} && date) &> >(tee {log})"

rule samtools_sort_index:
    input:
        os.path.join(mapping_dir("{assembly_type}"),
                     "{sample,[^\\.]+}.{read_tech,(hifi|ont)}.bam")
    output:
        bam = os.path.join(mapping_dir("{assembly_type}"),
                           "{sample,[^\\.]+}.{read_tech,(hifi|ont)}.sorted.bam"),
        bai = os.path.join(mapping_dir("{assembly_type}"),
                           "{sample,[^\\.]+}.{read_tech,(hifi|ont)}.sorted.bam.bai")
    threads:
        config["samtools"]["threads"]
    conda:
        os.path.join(ENV_DIR, "mapping.yaml")
    log:
        os.path.join(LOG_DIR, "samtools_sort/{assembly_type}/{sample}.{read_tech}.log")
    message:
        "samtools: sort+index {wildcards.assembly_type}/{wildcards.sample}.{wildcards.read_tech}"
    shell:
        "(date && samtools sort -@ {threads} -o {output.bam} {input} && "
        "samtools index -@ {threads} {output.bam} && date) &> >(tee {log})"

# --- Binning (MetaBAT2) ---
# NOTE: This uses -a with a space-separated BAM list as in your original draft.
# Standard MetaBAT2 expects a depth file (from jgi_summarize_bam_contig_depths).
# Keep as-is to mirror your prior behaviour.

def get_bams_for_binning(wc):
    assembly_type = wc.assembly_type
    bams = []
    if assembly_type == "hybrid_b26":
        # HiFi for the hybrid individual sample
        bams.append(os.path.join(mapping_dir(assembly_type),
                                 f"{HYBRID_INDIVIDUAL}.hifi.sorted.bam"))
        # ONT mapped BAM (sample name 'ont')
        bams.append(os.path.join(mapping_dir(assembly_type), "ont.ont.sorted.bam"))
    elif assembly_type == "hybrid_coassembly_bulk":
        for s in BULK_SAMPLES:
            bams.append(os.path.join(mapping_dir(assembly_type), f"{s}.hifi.sorted.bam"))
        bams.append(os.path.join(mapping_dir(assembly_type), "ont.ont.sorted.bam"))
    return bams

rule metabat2_binning:
    """Run MetaBAT2 on specified assemblies."""
    input:
        assembly = os.path.join(RESULTS_DIR, "assemblies/{assembly_type}/contigs.fasta"),
        bams = get_bams_for_binning
    output:
        out_prefix = os.path.join(binning_dir("{assembly_type}"), "{assembly_type}_bin"),
        done = touch(os.path.join(binning_dir("{assembly_type}"), "bins.DONE"))
    params:
        min_contig = config["metabat2"]["min_contig_length"],
        opts = config["metabat2"]["opts"],
        # NOTE: accept the *named* 'input' argument and join its 'bams' list
        bam_list = lambda wildcards, input: " ".join(input.bams)
    threads:
        config["metamdbg"]["threads"]
    conda:
        os.path.join(ENV_DIR, "metabat2.yaml")
    log:
        os.path.join(LOG_DIR, "metabat2/{assembly_type}.log")
    benchmark:
        os.path.join(BENCHMARK_DIR, "metabat2/{assembly_type}.txt")
    message:
        "MetaBAT2: binning {wildcards.assembly_type}"
    shell:
        "(date && metabat2 -i {input.assembly} -a {params.bam_list} "
        "-o {output.out_prefix} -m {params.min_contig} -t {threads} {params.opts} && date) &> >(tee {log})"

