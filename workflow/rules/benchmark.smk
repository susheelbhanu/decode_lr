##################################################
# BENCHMARK: ORF-caller + marker-gene comparison on a contig subset
#
# Questions:
#   1. prodigal-gv  vs  pyrodigal-gv           (ORF calling)
#   2. whole fasta  vs  20-way Split_Fasta     ("do chunks go faster?")
#   3. rpsblast/hmmsearch  vs  fetchMGs        (SCG / marker-gene calling)
#
# Existing commands (Split_Fasta, prodigal-gv, cat_orfs, hmmsearch/rpsblast,
# Filter_*, Extract_SCG) are COPIED VERBATIM from rules/annotation.smk; only
# the I/O wiring is adjusted. annotation.smk is left untouched.
# NEW rules: subset_contigs, pyrodigal-gv (whole + chunk), fetchMGs, summary.
#
# EXECUTION: submitted in parallel via profiles/slurm_annotation.
# Every rule asks for the SAME resources (uniform --cpus-per-task / --mem /
# partition) so the timing/memory numbers are comparable across jobs:
#   threads          = BENCH_THREADS   (uniform --cpus-per-task)
#   resources.mem_mb = BENCH_MEM_MB    (uniform --mem)
#   resources.slurm_partition = BENCH_PARTITION (ei-long; overrides the
#                      profile's ei-medium,ei-long,ei-largemem default so the
#                      ei-largemem submission trap is avoided)
# NOTE: qos still comes from the profile's cluster-config (slurm.yaml). Set
#       its __default__ qos to 'normal' (NOT qos-batch) or all jobs fail.
#
# Entry point: workflow/Snakefile_benchmark
##################################################

import os
import glob
from os.path import join, dirname, basename

# ------------------------------------------------------------------ #
#  Config  (reuses config/lr_config.yaml; override any with --config)
# ------------------------------------------------------------------ #
SAMPLE_KEY      = config.get("bench_sample", "hifi_individual_b10_d3_0200_con")
SUBSET_FRAC     = float(config.get("bench_fraction", 0.1))                 # 1/10th of contigs
SUBSET_SEED     = int(config.get("bench_seed", 11))
REPEATS         = int(config.get("bench_repeats", 3))                      # repeats on headline rules
SPLIT           = int(config["annotation"].get("prodigal_splits", 20))     # 20-way split, as pipeline

# ---- uniform resource request for EVERY job ----
BENCH_THREADS   = int(config.get("bench_threads", 8))      # uniform --cpus-per-task
BENCH_MEM_MB    = int(config.get("bench_mem_mb", 32000))   # uniform --mem
BENCH_PARTITION = config.get("bench_partition", "ei-long") # uniform partition (NOT ei-largemem)

OUT = os.path.abspath(config.get(
    "bench_dir",
    "/ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7/"
    "data/CEH_soil_project/lr_results_20260202/benchmark_orfcaller"))

# Containers -------------------------------------------------------- #
SIF_DIR        = "/ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7/data/CEH_soil_project/decode_lr/singularity_cache"
SIF_PRODIGALGV = join(SIF_DIR, "prodigal-gv_2.11.0.sif")   # existing
SIF_PYRODIGAL  = join(SIF_DIR, "pyrodigalgv_0.3.2.sif")    # NEW - pull first
SIF_FETCHMGS   = join(SIF_DIR, "fetchmgs.sif")             # NEW - pull first
SIF_HMM        = join(SIF_DIR, "gtdbtk_1.4.0.sif")         # existing (hmmsearch)
SIF_BLAST      = join(SIF_DIR, "blast_2.16.0.sif")         # existing (rpsblast)
SIF_PYENV      = join(SIF_DIR, "pythonenv_3.9.sif")        # existing (Split_Fasta, Filter_*, Extract_SCG)
SIF_SEQKIT     = join(SIF_DIR, "seqkit_2.6.1.sif")         # existing (subset)

# LongFlow scripts + SCG data  (same source as annotation.smk) ------ #
LONGFLOW_DIR = config["longflow_dir"]
SCRIPTS      = join(LONGFLOW_DIR, "snakenest", "python")
SCG_DATA     = join(LONGFLOW_DIR, "scg_data")

ANNOTATION_CFG = config["annotation"]

# SCG route: mirror annotation.smk exactly (cog_db -> rpsblast, else CheckM HMM)
USE_COG = "cog_db" in ANNOTATION_CFG
if USE_COG:
    COG_DB = ANNOTATION_CFG["cog_db"]
else:
    CHECKM_DB = ANNOTATION_CFG["checkm"]
    SCG_HMM   = f"{CHECKM_DB}/hmms/checkm.hmm"

SRC_GZ = config["assembly_paths"][SAMPLE_KEY]
TOOLS  = ["prodigal_gv", "pyrodigal_gv"]

wildcard_constraints:
    tool = "|".join(TOOLS),
    nb   = r"\d+",

# ================================================================== #
#  TARGET
# ================================================================== #
rule all:
    default_target: True
    input:
        join(OUT, "summary", "benchmark_summary.tsv"),
        join(OUT, "summary", "equivalence.tsv")
    output: touch("status/benchmark.done")

# ================================================================== #
#  NEW: subset 1/10 of the contigs  (shared input for everything)
# ================================================================== #
rule subset_contigs:
    input:  SRC_GZ
    output: join(OUT, "subset", "contigs.fa")
    params: frac=SUBSET_FRAC, seed=SUBSET_SEED
    threads: BENCH_THREADS
    resources:
        mem_mb          = BENCH_MEM_MB,
        slurm_partition = BENCH_PARTITION,
    benchmark: join(OUT, "benchmark", "subset_contigs.tsv")
    singularity: SIF_SEQKIT
    shell:
        "seqkit sample -p {params.frac} -s {params.seed} -o {output} {input}"

# ---- annotation.smk:`split_fasta` (verbatim command) ----
rule split_fasta:
    input:  join(OUT, "subset", "contigs.fa")
    output: expand(join(OUT, "split", "temp_splits", "Batch_{nb}"), nb=range(SPLIT))
    params: tmp = join(OUT, "split", "temp_splits")
    threads: BENCH_THREADS
    resources:
        mem_mb          = BENCH_MEM_MB,
        slurm_partition = BENCH_PARTITION,
    benchmark: join(OUT, "benchmark", "split_fasta.tsv")
    singularity: SIF_PYENV
    shell:
        "{SCRIPTS}/Split_Fasta.py {input} {SPLIT} -E -T {params.tmp}/Batch"

# ================================================================== #
#  ORF CALLING  -  WHOLE mode (single invocation on the full subset)
# ================================================================== #
# ARM 1a - prodigal-gv whole : command AS-IS from annotation.smk:`prodigal`
rule orf_prodigal_gv_whole:
    input:  join(OUT, "subset", "contigs.fa")
    output:
        faa = join(OUT, "prodigal_gv_whole", "contigs.faa"),
        fna = join(OUT, "prodigal_gv_whole", "contigs.fna"),
        gff = join(OUT, "prodigal_gv_whole", "contigs.gff"),
    threads: BENCH_THREADS                                # uniform request (tool is single-threaded)
    resources:
        mem_mb          = BENCH_MEM_MB,
        slurm_partition = BENCH_PARTITION,
    benchmark: repeat(join(OUT, "benchmark", "orf_prodigal_gv_whole.tsv"), REPEATS)
    singularity: SIF_PRODIGALGV
    shell:
        "prodigal-gv -i {input} -a {output.faa} -d {output.fna} -f gff -p meta -o {output.gff}"

# ARM 2a - pyrodigal-gv whole (NEW). -j uses the allocated cores (multithreaded)
rule orf_pyrodigal_gv_whole:
    input:  join(OUT, "subset", "contigs.fa")
    output:
        faa = join(OUT, "pyrodigal_gv_whole", "contigs.faa"),
        fna = join(OUT, "pyrodigal_gv_whole", "contigs.fna"),
        gff = join(OUT, "pyrodigal_gv_whole", "contigs.gff"),
    threads: BENCH_THREADS
    resources:
        mem_mb          = BENCH_MEM_MB,
        slurm_partition = BENCH_PARTITION,
    benchmark: repeat(join(OUT, "benchmark", "orf_pyrodigal_gv_whole.tsv"), REPEATS)
    singularity: SIF_PYRODIGAL
    shell:
        "pyrodigal-gv -i {input} -a {output.faa} -d {output.fna} -f gff -p meta -j {threads} -o {output.gff}"

# ================================================================== #
#  ORF CALLING  -  CHUNK mode (per Batch_{nb}, submitted in parallel, then cat)
# ================================================================== #
# ARM 1b - prodigal-gv per chunk : command AS-IS from annotation.smk:`prodigal`
rule orf_prodigal_gv_chunk:
    input:  join(OUT, "split", "temp_splits", "Batch_{nb}")
    output:
        faa = join(OUT, "prodigal_gv_chunk", "temp_splits", "Batch_{nb}.faa"),
        fna = join(OUT, "prodigal_gv_chunk", "temp_splits", "Batch_{nb}.fna"),
        gff = join(OUT, "prodigal_gv_chunk", "temp_splits", "Batch_{nb}.gff"),
    threads: BENCH_THREADS
    resources:
        mem_mb          = BENCH_MEM_MB,
        slurm_partition = BENCH_PARTITION,
    benchmark: join(OUT, "benchmark", "orf_prodigal_gv_chunk", "Batch_{nb}.tsv")
    singularity: SIF_PRODIGALGV
    shell:
        r"""
        if [ -s {input} ]; then
            prodigal-gv -i {input} -a {output.faa} -d {output.fna} -f gff -p meta -o {output.gff}
        else
            touch {output}
        fi
        """

# ARM 2b - pyrodigal-gv per chunk (NEW)
rule orf_pyrodigal_gv_chunk:
    input:  join(OUT, "split", "temp_splits", "Batch_{nb}")
    output:
        faa = join(OUT, "pyrodigal_gv_chunk", "temp_splits", "Batch_{nb}.faa"),
        fna = join(OUT, "pyrodigal_gv_chunk", "temp_splits", "Batch_{nb}.fna"),
        gff = join(OUT, "pyrodigal_gv_chunk", "temp_splits", "Batch_{nb}.gff"),
    threads: BENCH_THREADS
    resources:
        mem_mb          = BENCH_MEM_MB,
        slurm_partition = BENCH_PARTITION,
    benchmark: join(OUT, "benchmark", "orf_pyrodigal_gv_chunk", "Batch_{nb}.tsv")
    singularity: SIF_PYRODIGAL
    shell:
        r"""
        if [ -s {input} ]; then
            pyrodigal-gv -i {input} -a {output.faa} -d {output.fna} -f gff -p meta -j {threads} -o {output.gff}
        else
            touch {output}
        fi
        """

# ---- annotation.smk:`cat_orfs` (verbatim run block), one per caller ----
rule cat_prodigal_gv:
    input:
        faa = expand(join(OUT, "prodigal_gv_chunk", "temp_splits", "Batch_{nb}.faa"), nb=range(SPLIT)),
        fna = expand(join(OUT, "prodigal_gv_chunk", "temp_splits", "Batch_{nb}.fna"), nb=range(SPLIT)),
        gff = expand(join(OUT, "prodigal_gv_chunk", "temp_splits", "Batch_{nb}.gff"), nb=range(SPLIT)),
    output:
        faa = join(OUT, "prodigal_gv_chunk", "contigs.faa"),
        fna = join(OUT, "prodigal_gv_chunk", "contigs.fna"),
        gff = join(OUT, "prodigal_gv_chunk", "contigs.gff"),
    threads: BENCH_THREADS
    resources:
        mem_mb          = BENCH_MEM_MB,
        slurm_partition = BENCH_PARTITION,
    benchmark: join(OUT, "benchmark", "cat_prodigal_gv.tsv")
    run:
        shell("cat {input.faa} > {output.faa}")
        shell("cat {input.fna} > {output.fna}")
        header = next(open(input["gff"][0]))
        with open(output["gff"], "w") as handle_w:
            handle_w.write(header)
            for file in input["gff"]:
                if os.stat(file).st_size > 0:
                    with open(file) as handle_r:
                        _ = next(handle_r)
                        handle_w.writelines(line for line in handle_r)

rule cat_pyrodigal_gv:
    input:
        faa = expand(join(OUT, "pyrodigal_gv_chunk", "temp_splits", "Batch_{nb}.faa"), nb=range(SPLIT)),
        fna = expand(join(OUT, "pyrodigal_gv_chunk", "temp_splits", "Batch_{nb}.fna"), nb=range(SPLIT)),
        gff = expand(join(OUT, "pyrodigal_gv_chunk", "temp_splits", "Batch_{nb}.gff"), nb=range(SPLIT)),
    output:
        faa = join(OUT, "pyrodigal_gv_chunk", "contigs.faa"),
        fna = join(OUT, "pyrodigal_gv_chunk", "contigs.fna"),
        gff = join(OUT, "pyrodigal_gv_chunk", "contigs.gff"),
    threads: BENCH_THREADS
    resources:
        mem_mb          = BENCH_MEM_MB,
        slurm_partition = BENCH_PARTITION,
    benchmark: join(OUT, "benchmark", "cat_pyrodigal_gv.tsv")
    run:
        shell("cat {input.faa} > {output.faa}")
        shell("cat {input.fna} > {output.fna}")
        header = next(open(input["gff"][0]))
        with open(output["gff"], "w") as handle_w:
            handle_w.write(header)
            for file in input["gff"]:
                if os.stat(file).st_size > 0:
                    with open(file) as handle_r:
                        _ = next(handle_r)
                        handle_w.writelines(line for line in handle_r)

# ================================================================== #
#  SCG ANNOTATION  (existing route: rpsblast or hmmsearch)
#  Once per caller, off the WHOLE output. Commands AS-IS from annotation.smk.
# ================================================================== #
if USE_COG:
    rule scg_search:                                  # annotation.smk:`Batch_rpsblast`
        input:  join(OUT, "{tool}_whole", "contigs.faa")
        output: join(OUT, "{tool}_whole", "contigs.cogs.tsv")
        log:    join(OUT, "{tool}_whole", "contigs_cog.log")
        params: db = COG_DB
        threads: BENCH_THREADS
        resources:
            mem_mb          = BENCH_MEM_MB,
            slurm_partition = BENCH_PARTITION,
        benchmark: repeat(join(OUT, "benchmark", "scg_search_{tool}.tsv"), REPEATS)
        singularity: SIF_BLAST
        shell:
            "rpsblast -outfmt '6 qseqid sseqid evalue pident length slen qlen' "
            "-evalue 0.00001 -query {input} -db {params.db} -out {output} &>{log}"

    rule scg_parse:                                   # annotation.smk:`parse_cogs_annotation`
        input:  join(OUT, "{tool}_whole", "contigs.cogs.tsv")
        output: join(OUT, "{tool}_whole", "contigs_cogs_best_hits.tsv")
        threads: BENCH_THREADS
        resources:
            mem_mb          = BENCH_MEM_MB,
            slurm_partition = BENCH_PARTITION,
        benchmark: join(OUT, "benchmark", "scg_parse_{tool}.tsv")
        singularity: SIF_PYENV
        shell:
            "{SCRIPTS}/Filter_Cogs.py {input} --cdd_cog_file {SCG_DATA}/cdd_to_cog.tsv > {output}"
else:
    rule scg_search:                                  # annotation.smk:`hmmsearch`
        input:
            faa = join(OUT, "{tool}_whole", "contigs.faa"),
            db  = SCG_HMM,
        output: join(OUT, "{tool}_whole", "contigs_hmm.out")
        log:    join(OUT, "{tool}_whole", "contigs_hmm.log")
        threads: BENCH_THREADS
        resources:
            mem_mb          = BENCH_MEM_MB,
            slurm_partition = BENCH_PARTITION,
        benchmark: repeat(join(OUT, "benchmark", "scg_search_{tool}.tsv"), REPEATS)
        singularity: SIF_HMM
        shell:                                        # --cpu 1 kept verbatim from annotation.smk
            "hmmsearch --cut_tc --cpu 1 -o /dev/null --noali "
            "--domtblout {output} {input.db} {input.faa} >{log} 2>&1"

    rule scg_parse:                                   # annotation.smk:`parse_cogs_annotation`
        input:  join(OUT, "{tool}_whole", "contigs_hmm.out")
        output: join(OUT, "{tool}_whole", "contigs_cogs_best_hits.tsv")
        threads: BENCH_THREADS
        resources:
            mem_mb          = BENCH_MEM_MB,
            slurm_partition = BENCH_PARTITION,
        benchmark: join(OUT, "benchmark", "scg_parse_{tool}.tsv")
        singularity: SIF_PYENV
        shell:
            "{SCRIPTS}/Filter_scg_hmm.py {input} --cog_hmm {SCG_DATA}/scg_hmm_selected.txt {output}"

# ---- annotation.smk:`extract_SCG_sequences` (verbatim command) ----
rule extract_scg:
    input:
        annotation = join(OUT, "{tool}_whole", "contigs_cogs_best_hits.tsv"),
        gff        = join(OUT, "{tool}_whole", "contigs.gff"),
        fna        = join(OUT, "{tool}_whole", "contigs.fna"),
    output: join(OUT, "{tool}_whole", "contigs_SCG.fna")
    threads: BENCH_THREADS
    resources:
        mem_mb          = BENCH_MEM_MB,
        slurm_partition = BENCH_PARTITION,
    benchmark: join(OUT, "benchmark", "extract_scg_{tool}.tsv")
    singularity: SIF_PYENV
    shell:
        "{SCRIPTS}/Extract_SCG.py {input.fna} {input.annotation} "
        "{SCG_DATA}/scg_cogs_min0.97_max1.03_unique_genera.txt {input.gff} > {output}"

# ================================================================== #
#  NEW: fetchMGs marker-gene extraction (alternative SCG route)
#  Runs off the pyrodigal-gv proteins+genes (mOTU 40 universal MGs).
#  NOTE: '-x ""' = find hmmsearch/seqtk on $PATH (present in the container).
#  If it can't locate its bundled HMM library, add '-l <fetchMGs lib path>'.
#  Confirm entrypoint: singularity exec <sif> which fetchMGs.pl
# ================================================================== #
rule fetchmg_pyrodigal_gv:
    input:
        faa = join(OUT, "pyrodigal_gv_whole", "contigs.faa"),
        fna = join(OUT, "pyrodigal_gv_whole", "contigs.fna"),
    output:
        done = join(OUT, "pyrodigal_gv_whole", "fetchMGs.done"),
    params:
        outdir = join(OUT, "pyrodigal_gv_whole", "fetchMGs_out"),
    threads: BENCH_THREADS
    resources:
        mem_mb          = BENCH_MEM_MB,
        slurm_partition = BENCH_PARTITION,
    benchmark: repeat(join(OUT, "benchmark", "fetchmg_pyrodigal_gv.tsv"), REPEATS)
    singularity: SIF_FETCHMGS
    shell:
        r"""
        rm -rf {params.outdir} && mkdir -p {params.outdir}
        fetchMGs extraction {input.faa} gene {params.outdir} -d {input.fna} -t {threads}
        touch {output.done}
        """

# ================================================================== #
#  SUMMARY  (collate benchmark TSVs; ORF/SCG/MG counts)   [run locally]
# ================================================================== #
rule summary:
    input:
        whole = expand(join(OUT, "{tool}_whole", "contigs.faa"), tool=TOOLS),
        chunk = expand(join(OUT, "{tool}_chunk", "contigs.faa"), tool=TOOLS),
        scg   = expand(join(OUT, "{tool}_whole", "contigs_SCG.fna"), tool=TOOLS),
        mg    = join(OUT, "pyrodigal_gv_whole", "fetchMGs.done"),
    output:
        bench = join(OUT, "summary", "benchmark_summary.tsv"),
        equiv = join(OUT, "summary", "equivalence.tsv"),
    params:
        bdir   = join(OUT, "benchmark"),
        outdir = OUT,
        tools  = TOOLS,
        mgdir  = join(OUT, "pyrodigal_gv_whole", "fetchMGs_out"),
    run:
        os.makedirs(dirname(output.bench), exist_ok=True)

        # ---- collate every benchmark/**/*.tsv (repeat() -> one row per run) ----
        with open(output.bench, "w") as out:
            out.write("step\trun\ts_seconds\th_m_s\tmax_rss_MB\tmax_pss_MB\t"
                      "mean_load\tcpu_time\n")
            for f in sorted(glob.glob(join(params.bdir, "**", "*.tsv"), recursive=True)):
                step = os.path.relpath(f, params.bdir)[:-4]   # keeps chunk subdir in the name
                with open(f) as fh:
                    header = fh.readline().rstrip("\n").split("\t")
                    run_i = 0
                    for line in fh:
                        vals = line.rstrip("\n").split("\t")
                        if not vals or vals == [""]:
                            continue
                        run_i += 1
                        d = dict(zip(header, vals))
                        out.write("\t".join([
                            step, str(run_i),
                            d.get("s", ""), d.get("h:m:s", ""),
                            d.get("max_rss", ""), d.get("max_pss", ""),
                            d.get("mean_load", ""), d.get("cpu_time", ""),
                        ]) + "\n")

        # ---- equivalence: ORF + SCG + fetchMGs counts ----
        def count_fasta(p):
            n = 0
            if os.path.exists(p):
                with open(p) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            n += 1
            return n

        with open(output.equiv, "w") as out:
            out.write("method\tn_orfs\tn_scg_or_mg_seqs\n")
            for t in params.tools:
                faa = join(params.outdir, t + "_whole", "contigs.faa")
                scg = join(params.outdir, t + "_whole", "contigs_SCG.fna")
                out.write(f"{t}\t{count_fasta(faa)}\t{count_fasta(scg)}\n")
            mg_faa = sorted(glob.glob(join(params.mgdir, "*.faa")))
            mg_n = sum(count_fasta(p) for p in mg_faa)
            pyr = join(params.outdir, "pyrodigal_gv_whole", "contigs.faa")
            out.write(f"pyrodigal_gv+fetchMGs\t{count_fasta(pyr)}\t{mg_n}\n")

