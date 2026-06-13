##################################################
# ANNOTATION UPDATED: LongFlow annotation workflow for
# pre-computed decode_lr assemblies
#
# Changes from annotation.smk (all original rules kept as-is):
#   1. pyrodigal       — whole-fasta ORF calling, replaces split→prodigal→cat chain
#   2. kofamscan_direct — KO annotation direct on contigs.faa (step 2 priority)
#   3. parse_scg_whole — single-file SCG parse, replaces batch-based parse_cogs_annotation
#   4. fetchMGs        — parallel SCG route (mOTU 40 universal marker genes)
#   Ruleorder directives resolve ambiguity with original batch-based rules.
#   CAT_annotation still uses the original split_fasta→prodigal batch chain.
#
# New config keys needed under annotation:
#   sif_pyrodigal: /hpc-home/kar23heg/singularity_cache/pyrodigal-gv.sif
#   sif_fetchmgs:  /hpc-home/kar23heg/singularity_cache/fetchmgs.sif
#
# Entry point: workflow/Snakefile_annotation
##################################################

import os
import glob
import re
from functools import partial
from collections import defaultdict, Counter
from os.path import basename, dirname, realpath, abspath
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp

# ------------------------------------------------------------------ #
#  SLURM resource setup (mirrors init.smk without the sample loading)
# ------------------------------------------------------------------ #

def _build_slurm_partitions(cfg):
    parts = []
    try:
        for _, specs in cfg.get("slurm_partitions", {}).items():
            name       = specs.get("name", "")
            min_mem_mb = 1000 * int(specs.get("min_mem", 0))
            max_mem_mb = 1000 * int(specs.get("max_mem", 0))
            min_thr    = int(specs.get("min_threads", 0))
            max_thr    = int(specs.get("max_threads", 0))
            parts.append([name, min_mem_mb, max_mem_mb, min_thr, max_thr])
    except Exception:
        parts = []
    return parts

SLURM_PARTITIONS = _build_slurm_partitions(config)

def get_resource_real(wildcards, input, threads, attempt,
                      SLURM_PARTITIONS="", mode="", mult=2, min_size=10000):
    def _return(mem, partition, thr, mode):
        if mode == "mem":       return int(mem)
        if mode == "partition": return partition
        if mode == "threads":   return int(thr)
        return int(mem)

    try:
        mem = max((input.size // 1_000_000) * attempt * mult,
                  attempt * min_size * mult)
    except Exception:
        mem = attempt * min_size * mult

    if not SLURM_PARTITIONS or SLURM_PARTITIONS[0][0] == "":
        return _return(mem, "", threads, mode)

    mem_ok = [p for p in SLURM_PARTITIONS if p[2] >= mem] \
             or [max(SLURM_PARTITIONS, key=lambda x: x[2])]
    thr_ok = [p for p in mem_ok if p[4] >= threads] \
             or [max(mem_ok, key=lambda x: x[4])]
    name, min_mem, max_mem, min_thr, max_thr = min(
        thr_ok, key=lambda x: [x[1], x[3], x[2], x[4]]
    )
    mem_final = min(max(mem, min_mem), max_mem)
    thr_final = min(max(threads, min_thr), max_thr)
    return _return(mem_final, name, thr_final, mode)

def get_resource(mode, **kwargs):
    return partial(get_resource_real, SLURM_PARTITIONS=SLURM_PARTITIONS,
                   mode=mode, **kwargs)

# ------------------------------------------------------------------ #
#  Paths to LongFlow helper scripts and SCG data
# ------------------------------------------------------------------ #
LONGFLOW_DIR = config["longflow_dir"]
SCRIPTS      = os.path.join(LONGFLOW_DIR, "snakenest", "python")
SCG_DATA     = os.path.join(LONGFLOW_DIR, "scg_data")

# Output directory for annotation results
ANNOT_DIR = os.path.abspath(config["annotation_dir"])

# ------------------------------------------------------------------ #
#  Annotation config
# ------------------------------------------------------------------ #
ANNOTATION_CFG  = config["annotation"]
PRODIGAL_SPLIT  = int(ANNOTATION_CFG.get("prodigal_splits", 20))
PID_ORF         = ANNOTATION_CFG.get("orf_dereplication", 0.99)

DIAMOND         = ANNOTATION_CFG.get("diamond", {})
BLCA            = ANNOTATION_CFG.get("blca", {"db": None, "taxa": None, "bin": None})
KO_HMM          = ANNOTATION_CFG.get("kofamscan", {}).get("profiles", "")
KO_HMM_CUTOFFS  = ANNOTATION_CFG.get("kofamscan", {}).get("ko_list", "")
SMK_CONFIG      = ANNOTATION_CFG.get("kofamscan", {}).get("smk_config", "")
GENOMAD_DB      = ANNOTATION_CFG.get("genomad_db", "")
CAT_DB          = ANNOTATION_CFG.get("cat_db", "")
CAT_PATH        = ANNOTATION_CFG.get("cat_path", "")
IP_DB           = ANNOTATION_CFG.get("ip_db", "")
CHECKM2_DB      = ANNOTATION_CFG.get("checkm2", "")

_SIF_CACHE    = "/ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7/data/CEH_soil_project/decode_lr/singularity_cache"
SIF_PYRODIGAL  = ANNOTATION_CFG.get("sif_pyrodigal",  f"{_SIF_CACHE}/pyrodigalgv_0.3.2.sif")
SIF_FETCHMGS   = ANNOTATION_CFG.get("sif_fetchmgs",   f"{_SIF_CACHE}/fetchmgs.sif")
SIF_PYTHONENV  = ANNOTATION_CFG.get("sif_pythonenv",  f"{_SIF_CACHE}/pythonenv_3.9.sif")
SIF_PRODIGAL   = ANNOTATION_CFG.get("sif_prodigal",   f"{_SIF_CACHE}/prodigal-gv_2.11.0.sif")
SIF_GTDBTK     = ANNOTATION_CFG.get("sif_gtdbtk",     f"{_SIF_CACHE}/gtdbtk_1.4.0.sif")
SIF_BARRNAP    = ANNOTATION_CFG.get("sif_barrnap",    f"{_SIF_CACHE}/barrnap_0.9.sif")
SIF_BLAST      = ANNOTATION_CFG.get("sif_blast",      f"{_SIF_CACHE}/blast_2.16.0.sif")
SIF_DIAMOND    = ANNOTATION_CFG.get("sif_diamond",    f"{_SIF_CACHE}/diamond_2.0.9.sif")
SIF_MMSEQS2    = ANNOTATION_CFG.get("sif_mmseqs2",    f"{_SIF_CACHE}/mmseqs2_latest.sif")
SIF_GENOMAD    = ANNOTATION_CFG.get("sif_genomad",    f"{_SIF_CACHE}/genomad_1.8.0.sif")

# SCG detection: use cog_db if present, else checkm HMM
if "cog_db" in ANNOTATION_CFG:
    COG_DB  = ANNOTATION_CFG["cog_db"]
    SCG     = {
        line.rstrip().split("\t")[0]
        for line in open(f"{SCG_DATA}/scg_cogs_min0.97_max1.03_unique_genera.txt")
    }
else:
    CHECKM_DB = ANNOTATION_CFG["checkm"]
    SCG_HMM   = f"{CHECKM_DB}/hmms/checkm.hmm"
    SCG       = {
        line.rstrip().split("\t")[0]
        for line in open(f"{SCG_DATA}/scg_hmm_selected.txt")
        if line.rstrip().split("\t")[2] == "fine"
    }

NB_SCG = float(len(SCG))

# ------------------------------------------------------------------ #
#  Assembly paths
# ------------------------------------------------------------------ #
ASSEMBLY_GZ  = config["assembly_paths"]
ASSEMBLIES   = list(ASSEMBLY_GZ.keys())

# ------------------------------------------------------------------ #
#  Build annotation output list
# ------------------------------------------------------------------ #
ANNOTATION_OUTPUT = [
    f"SCG_cluster_{PID_ORF}.tsv",
    f"summary_{PID_ORF}.tsv",
]
if DIAMOND:
    for annot in DIAMOND:
        ANNOTATION_OUTPUT.append(f"contigs_{annot}_best_hits.tsv")
if BLCA.get("db"):
    ANNOTATION_OUTPUT.append("16S_blca.tsv")
if KO_HMM:
    ANNOTATION_OUTPUT.append("contigs_KEGG_best_hits.tsv")
if GENOMAD_DB:
    ANNOTATION_OUTPUT.append("genomad/genomad.done")
if CAT_DB:
    ANNOTATION_OUTPUT.append("CAT_contigs_taxonomy.tsv")
if IP_DB:
    ANNOTATION_OUTPUT.append("contigs_IP_best_hits.tsv")

# fetchMGs always added when SIF is configured
if SIF_FETCHMGS:
    ANNOTATION_OUTPUT.append("fetchMGs.done")

# ------------------------------------------------------------------ #
#  Wildcard constraints
# ------------------------------------------------------------------ #
wildcard_constraints:
    assembly   = "|".join(re.escape(a) for a in ASSEMBLIES),
    annotation = "|".join(DIAMOND.keys()) if DIAMOND else "NO_ANNOTATION",

# ------------------------------------------------------------------ #
#  Rule priority: new whole-fasta rules take precedence over
#  the original batch-based rules where outputs overlap.
# ------------------------------------------------------------------ #

# ================================================================== #
#  TARGET RULE
# ================================================================== #

rule annotation_all:
    input:  expand("{annot_dir}/{assembly}/annotation/{output}", annot_dir=ANNOT_DIR, assembly=ASSEMBLIES, output=ANNOTATION_OUTPUT)
    output: touch(ANNOT_DIR + "/.annotation_all.done")

# ================================================================== #
#  PREPARE CONTIGS  (decompress existing decode_lr assembly)
# ================================================================== #

rule prepare_contigs:
    """Decompress the pre-computed assembly into the working tree."""
    input:  lambda w: ASSEMBLY_GZ[w.assembly]
    output: "{annot_dir}/{assembly}/contigs/contigs.fa"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb          = get_resource("mem"),
    shell: "gunzip -c {input} > {output}"

# ================================================================== #
#  STEP 1 — PYRODIGAL  (whole-fasta, replaces split→prodigal→cat)
# ================================================================== #

rule pyrodigal:
    """Whole-fasta ORF prediction with pyrodigal-gv (multithreaded)."""
    input:  "{annot_dir}/{assembly}/contigs/contigs.fa"
    output:
        faa = "{annot_dir}/{assembly}/annotation/contigs.faa",
        fna = "{annot_dir}/{assembly}/annotation/contigs.fna",
        gff = "{annot_dir}/{assembly}/annotation/contigs.gff",
    priority: 100
    threads: 32
    resources:
        slurm_partition = get_resource("partition", min_size=150000),
        mem_mb          = get_resource("mem", min_size=150000),
    singularity: SIF_PYRODIGAL
    shell: """
    pyrodigal-gv -i {input} -a {output.faa} -d {output.fna} \
        -f gff -p meta -j {threads} -o {output.gff}
    """

# ================================================================== #
#  ORIGINAL: split_fasta / prodigal / cat_orfs
#  Kept as-is — still used by CAT_annotation (needs nucleotide batches)
# ================================================================== #

rule split_fasta:
    input:  "{annot_dir}/{assembly}/contigs/contigs.fa"
    output: expand("{{annot_dir}}/{{assembly}}/annotation/temp_splits/Batch_{nb}", nb=range(PRODIGAL_SPLIT))
    params: tmp = "{annot_dir}/{assembly}/annotation/temp_splits"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb          = get_resource("mem"),
    singularity: SIF_PYTHONENV
    shell: "{SCRIPTS}/Split_Fasta.py {input} {PRODIGAL_SPLIT} -E -T {params.tmp}/Batch"

rule prodigal:
    input:  "{path}/temp_splits/Batch_{nb}"
    output:
        faa = "{path}/temp_splits/Batch_{nb}.faa",
        fna = "{path}/temp_splits/Batch_{nb}.fna",
        gff = "{path}/temp_splits/Batch_{nb}.gff",
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb          = get_resource("mem"),
    singularity: SIF_PRODIGAL
    shell: """
    if [ -s {input} ]; then
        prodigal-gv -i {input} -a {output.faa} -d {output.fna} -f gff -p meta -o {output.gff}
    else
        touch {output}
    fi
    """

# rule cat_orfs — DISABLED: pyrodigal produces contigs.faa/fna/gff directly
# rule cat_orfs:
#     input:
#         faa = expand("{{path}}/temp_splits/Batch_{nb}.faa", nb=range(PRODIGAL_SPLIT)),
#         fna = expand("{{path}}/temp_splits/Batch_{nb}.fna", nb=range(PRODIGAL_SPLIT)),
#         gff = expand("{{path}}/temp_splits/Batch_{nb}.gff", nb=range(PRODIGAL_SPLIT)),
#     output:
#         faa = "{path}/contigs.faa",
#         fna = "{path}/contigs.fna",
#         gff = "{path}/contigs.gff",
#     resources:
#         slurm_partition = get_resource("partition"),
#         mem_mb          = get_resource("mem"),
#     run:
#         shell("cat {input.faa} > {output.faa}")
#         shell("cat {input.fna} > {output.fna}")
#         header = next(open(input["gff"][0]))
#         with open(output["gff"], "w") as handle_w:
#             handle_w.write(header)
#             for file in input["gff"]:
#                 if os.stat(file).st_size > 0:
#                     with open(file) as handle_r:
#                         _ = next(handle_r)
#                         handle_w.writelines(line for line in handle_r)

# ================================================================== #
#  ORIGINAL: SCG ANNOTATION  (COG-db path or CheckM HMM path)
#  hmmsearch / Batch_rpsblast run on individual batch .faa files.
#  parse_cogs_annotation aggregates all batches — still used by CAT path.
# ================================================================== #

if "cog_db" in ANNOTATION_CFG:
    rule Batch_rpsblast:
        input:   "{path}.faa"
        output:  "{path}.cogs.tsv"
        params:  db = COG_DB
        log:     "{path}_cog.log"
        resources:
            slurm_partition = get_resource("partition"),
            mem_mb          = get_resource("mem"),
        singularity: SIF_BLAST
        shell: """
            rpsblast -outfmt '6 qseqid sseqid evalue pident length slen qlen' \
                -evalue 0.00001 -query {input} -db {params.db} -out {output} &>{log}
        """

    # rule parse_cogs_annotation — DISABLED: parse_scg_whole replaces the batch aggregation
    # rule parse_cogs_annotation:
    #     input:  Batch = expand("{{path}}/annotation/temp_splits/Batch_{nb}.cogs.tsv", nb=range(PRODIGAL_SPLIT))
    #     output:
    #         cog = "{path}/annotation/contigs_cogs_best_hits.tsv",
    #         cat = temp("{path}/annotation/contigs_Cog.out"),
    #     resources:
    #         slurm_partition = get_resource("partition"),
    #         mem_mb          = get_resource("mem"),
    #     singularity: SIF_PYTHONENV
    #     shell: """
    #         cat {input} > {output.cat}
    #         {SCRIPTS}/Filter_Cogs.py {output.cat} \
    #             --cdd_cog_file {SCG_DATA}/cdd_to_cog.tsv > {output.cog}
    #     """

else:
    rule hmmsearch:
        input:
            faa = "{path}.faa",
            db  = SCG_HMM,
        output: "{path}_hmm.out"
        log:    "{path}_hmm.log"
        resources:
            slurm_partition = get_resource("partition", min_size=150000),
            mem_mb          = get_resource("mem", min_size=150000),
        singularity: SIF_GTDBTK
        shell: """
        if [ ! -s {input.faa} ]; then
            touch {output}
        else
            hmmsearch --cut_tc --cpu 1 -o /dev/null --noali \
                --domtblout {output} {input.db} {input.faa} >{log} 2>&1
        fi
        """

    # rule parse_cogs_annotation — DISABLED: parse_scg_whole replaces the batch aggregation
    # rule parse_cogs_annotation:
    #     input:  Batch = expand("{{path}}/temp_splits/Batch_{nb}_hmm.out", nb=range(PRODIGAL_SPLIT))
    #     output:
    #         cat = temp("{path}/contigs_hmm.out"),
    #         cog = "{path}/contigs_cogs_best_hits.tsv",
    #     resources:
    #         slurm_partition = get_resource("partition"),
    #         mem_mb          = get_resource("mem"),
    #     singularity: SIF_PYTHONENV
    #     shell: """
    #         cat {input} > {output.cat}
    #         {SCRIPTS}/Filter_scg_hmm.py {output.cat} \
    #             --cog_hmm {SCG_DATA}/scg_hmm_selected.txt {output.cog}
    #     """

# ================================================================== #
#  STEP 3 — PARSE SCG WHOLE  (single-file, for pyrodigal whole path)
#  hmmsearch on contigs.faa produces contigs_hmm.out directly;
#  this rule parses it to contigs_cogs_best_hits.tsv.
#  ruleorder above ensures this fires instead of parse_cogs_annotation.
# ================================================================== #

if "cog_db" in ANNOTATION_CFG:
    rule parse_scg_whole:
        input: "{path}/annotation/contigs.cogs.tsv"
        output: "{path}/annotation/contigs_cogs_best_hits.tsv"
        priority: 80
        resources:
            slurm_partition = get_resource("partition"),
            mem_mb          = get_resource("mem"),
        singularity: SIF_PYTHONENV
        shell: """
            {SCRIPTS}/Filter_Cogs.py {input} \
                --cdd_cog_file {SCG_DATA}/cdd_to_cog.tsv > {output}
        """
else:
    rule parse_scg_whole:
        input: "{path}/annotation/contigs_hmm.out"
        output: "{path}/annotation/contigs_cogs_best_hits.tsv"
        priority: 80
        resources:
            slurm_partition = get_resource("partition"),
            mem_mb          = get_resource("mem"),
        singularity: SIF_PYTHONENV
        shell: """
            {SCRIPTS}/Filter_scg_hmm.py {input} \
                --cog_hmm {SCG_DATA}/scg_hmm_selected.txt {output}
        """

# ================================================================== #
#  ORIGINAL: SCG EXTRACTION AND CLUSTERING  (unchanged, patterns match)
# ================================================================== #

rule extract_SCG_sequences:
    input:
        annotation = "{filename}_cogs_best_hits.tsv",
        gff        = "{filename}.gff",
        fna        = "{filename}.fna",
    output: "{filename}_SCG.fna"
    resources:
        slurm_partition = get_resource("partition", min_size=150000),
        mem_mb          = get_resource("mem", min_size=150000),
    shell: """
        {SCRIPTS}/Extract_SCG.py {input.fna} {input.annotation} \
            {SCG_DATA}/scg_cogs_min0.97_max1.03_unique_genera.txt \
            {input.gff} > {output}
    """

rule cluster_SCG:
    input:  "{path}/contigs_SCG.fna"
    output: "{path}/contigs_{pid}_mmseqs_cluster.tsv"
    params:
        tmp = "{path}/mmseq_tmp_{pid}",
        out = "{path}/contigs_{pid}_mmseqs",
    threads: 10
    wildcard_constraints:
        pid = r"0\.\d+|1\.0+"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb          = get_resource("mem"),
    singularity: SIF_MMSEQS2
    shell: """
        mmseqs easy-cluster {input} {params.out} {params.tmp} \
            --min-seq-id {wildcards.pid} -c 0.8 \
            --cov-mode 1 --alignment-mode 3 --threads {threads}
    """

rule scg_cluster_def:
    input:
        clu = "{path}/contigs_{pid}_mmseqs_cluster.tsv",
        scg = "{path}/contigs_SCG.fna",
    output: "{path}/SCG_cluster_{pid}.tsv"
    run:
        orf_to_scg  = {header.split()[0]: header.split()[1] for header, _ in sfp(open(input["scg"]))}
        orf_to_clu  = defaultdict(lambda: defaultdict(list))
        for line in open(input["clu"]):
            rep, orf = line.rstrip().split("\t")
            orf_to_clu[orf_to_scg[rep]][rep].append(orf)
        with open(output[0], "w") as handle:
            handle.writelines(
                "%s\t%s\t%s\n" % (cog, index, "\t".join(orfs))
                for cog, clus in orf_to_clu.items()
                for index, (ref, orfs) in enumerate(clus.items())
            )

rule cluster_ORF:
    input:  "{path}/contigs.fna"
    output: "{path}/contigs_orfs_{pid}_mmseqs_cluster.tsv"
    params:
        tmp = "{path}/mmseq_tmp_orfs_{pid}",
        out = "{path}/contigs_orfs_{pid}_mmseqs",
    threads: 10
    resources:
        slurm_partition = get_resource("partition", mult=8),
        mem_mb          = get_resource("mem", mult=8),
    singularity: SIF_MMSEQS2
    shell: """
        mmseqs easy-cluster {input} {params.out} {params.tmp} \
            --min-seq-id {wildcards.pid} -c 0.8 \
            --cov-mode 1 --alignment-mode 3 --threads {threads}
    """

# ================================================================== #
#  ORIGINAL: CONTIG QUALITY  (unchanged)
# ================================================================== #

rule bogus_bed:
    input:  "{path}/contigs/contigs.fa"
    output: temp("{path}/annotation/contigs.bedtemp")
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb          = get_resource("mem"),
    shell: "{SCRIPTS}/bogus_bed.py -i {input} -o {output}"

rule sort_bed:
    input:
        bed  = "{path}/annotation/{type}.bedtemp",
        cont = "{path}/contigs/contigs.fa",
    output:
        bed   = "{path}/annotation/{type}.bed",
        gfile = "{path}/annotation/{type}_bedtools_target_definition.tsv",
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb          = get_resource("mem"),
    shell: "{SCRIPTS}/sort_bed.py {input.bed} {input.cont} {output.bed} -g {output.gfile}"

rule get_component_quality:
    input:
        scgs     = "{path}/annotation/SCG_cluster_{pid}.tsv",
        cogs     = "{path}/annotation/contigs_cogs_best_hits.tsv",
        cont_len = "{path}/annotation/contigs.bed",
    output: summary = "{path}/annotation/summary_{pid}.tsv"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb          = get_resource("mem"),
    run:
        contig_to_len = {
            line.rstrip().split("\t")[0]: int(line.rstrip().split("\t")[2])
            for line in open(input["cont_len"])
        }
        orfs_to_cog  = {}
        cog_to_orfs  = defaultdict(list)
        for line in open(input["scgs"]):
            cog = "_".join(line.rstrip().split("\t")[0:2])
            for orf in line.rstrip().split("\t")[2:]:
                orfs_to_cog[orf] = cog
                cog_to_orfs[cog].append(orf)
        scg_to_cog = {
            line.rstrip().split("\t")[0]: line.rstrip().split("\t")[1]
            for index, line in enumerate(open(input["cogs"])) if index > 0
        }
        contigs_to_scgs = defaultdict(set)
        for orf, cog in orfs_to_cog.items():
            contig = "_".join(orf.split("_")[:-1])
            contigs_to_scgs[contig].add(cog)

        get_cogs      = lambda x: Counter([el.split("_")[0] for el in contigs_to_scgs[x]])
        completion    = lambda x: sum(v >= 1 for v in get_cogs(x).values()) / NB_SCG
        contamination = lambda x: (sum(get_cogs(x).values()) - len(get_cogs(x))) / NB_SCG

        sorted_contigs = sorted(contig_to_len.keys(), key=lambda x: -contig_to_len[x])
        with open(output["summary"], "w") as handle:
            handle.write("contig\tcomp\tcont\tnuc_size\n")
            handle.writelines(
                "%s\t%s\t%s\t%s\n" % (
                    contig,
                    "{:.0%}".format(completion(contig)),
                    "{:.0%}".format(contamination(contig)),
                    contig_to_len[contig],
                )
                for contig in sorted_contigs
            )

# ================================================================== #
#  ORIGINAL: DIAMOND ANNOTATION  (unchanged)
# ================================================================== #

rule diamond:
    input:   "{filename}.faa"
    output:  "{filename}_{annotation}.m8"
    params:  db  = lambda w: DIAMOND[w.annotation]["db"]
    log:     "{filename}_{annotation}_diamond.log"
    threads: 32
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb          = get_resource("mem"),
    singularity: SIF_DIAMOND
    shell: """
        diamond blastp --more-sensitive \
            -d {params.db} -q {input} -p {threads} -o {output} \
            -f6 qseqid sseqid qstart qend qlen sstart send slen length pident evalue bitscore \
            &>{log}
    """

rule annotation_diamond:
    input:   "{path}/annotation/{filename}_{annotation}.m8"
    output:  "{path}/annotation/{filename}_{annotation}_best_hits.tsv"
    params:
        annotation  = lambda w: DIAMOND[w.annotation]["annotation"],
        Bitscore    = lambda w: DIAMOND[w.annotation]["filter"][0],
        Evalue      = lambda w: DIAMOND[w.annotation]["filter"][1],
        PID         = lambda w: DIAMOND[w.annotation]["filter"][2],
        subject_pid = lambda w: DIAMOND[w.annotation]["filter"][3],
        subject_cov = lambda w: DIAMOND[w.annotation]["filter"][4],
        query_cov   = lambda w: DIAMOND[w.annotation]["filter"][5],
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb          = get_resource("mem"),
    shell: """
        {SCRIPTS}/M8_Filtering.py {input} \
            -D {params.annotation} \
            -B {params.Bitscore} -E {params.Evalue} -P {params.PID} \
            -R {params.subject_pid} -C {params.subject_cov} -Q {params.query_cov} \
            > {output}
    """

# ================================================================== #
#  STEP 2 — KO ANNOTATION  (hmmsearch against KO profiles, whole faa)
#  Mirrors ko_annotation.snake logic but on whole contigs.faa.
#  Two rules: hmmsearch (cluster job) → parse (local Python run:).
# ================================================================== #

# Parse per-KO bitscore/domain cutoffs from the ko_list file
def _parse_ko_cutoffs(ko_list_path):
    cutoffs = {}
    with open(ko_list_path) as fh:
        _ = fh.readline()
        for line in fh:
            parts = line.strip().split("\t")
            if parts[1] == "-":
                cutoffs[parts[0]] = [0.0, "NA"]
            else:
                cutoffs[parts[0]] = [float(parts[1]), parts[2]]
    return cutoffs

def _parse_ko_domtblout(path, ko_cutoffs, evalue_cutoff=1e-5):
    from collections import defaultdict
    hits = defaultdict(list)
    for line in open(path):
        if line.startswith("#"):
            continue
        cols = line.strip().split()[:22]
        if len(cols) < 22:
            continue
        seqname, _, tlen, ko_hmm, _, qlen, _, score, \
        _, _, _, _, i_evalue, dom_score, _, _, \
        _, seq_from, seq_to, _, _, _ = cols
        tlen, qlen = int(tlen), int(qlen)
        span = range(int(seq_from), int(seq_to))
        if ko_hmm in ko_cutoffs:
            mode = ko_cutoffs[ko_hmm][1]
            threshold = ko_cutoffs[ko_hmm][0]
            if mode == "full" and float(score) < threshold:
                continue
            elif mode == "domain" and float(dom_score) < threshold:
                continue
        if float(i_evalue) > evalue_cutoff:
            continue
        hits[seqname].append((ko_hmm, float(i_evalue), span))
    # keep best non-overlapping hit per orf
    result = {}
    for orf, annots in hits.items():
        if len(annots) == 1:
            result[orf] = [annots[0][0]]
            continue
        annots_sorted = sorted(annots, key=lambda x: x[1])
        chosen, excluded = [], set()
        for ko, ev, rng in annots_sorted:
            if ko in excluded:
                continue
            chosen.append(ko)
            for ko2, ev2, rng2 in annots_sorted:
                if set(rng).intersection(rng2):
                    excluded.add(ko2)
        result[orf] = chosen
    return result

if KO_HMM:
    rule kofamscan_direct:
        """hmmsearch against KO HMM profiles on whole contigs.faa."""
        input:  faa = "{path}/annotation/contigs.faa"
        output: domtbl = "{path}/annotation/contigs_ko.out"
        priority: 90
        threads: 32
        resources:
            slurm_partition = get_resource("partition", min_size=150000),
            mem_mb          = get_resource("mem", min_size=150000),
        singularity: SIF_GTDBTK
        shell: """
        if [ -s {input.faa} ]; then
            hmmsearch --cpu {threads} -o /dev/null --noali \
                --domtblout {output.domtbl} {KO_HMM} {input.faa}
        else
            touch {output.domtbl}
        fi
        """

    localrules: kofamscan_parse
    rule kofamscan_parse:
        """Parse KO domtblout → contigs_KEGG_best_hits.tsv."""
        input:  domtbl = "{path}/annotation/contigs_ko.out"
        output: tsv    = "{path}/annotation/contigs_KEGG_best_hits.tsv"
        priority: 90
        run:
            ko_cutoffs = _parse_ko_cutoffs(KO_HMM_CUTOFFS)
            ko_to_def  = {
                line.rstrip().split("\t")[0]: line.rstrip().split("\t")[-1]
                for line in open(KO_HMM_CUTOFFS)
                if not line.startswith("knum")
            }
            orf_to_kos = _parse_ko_domtblout(input.domtbl, ko_cutoffs)
            with open(output.tsv, "w") as fh:
                fh.write("orf\tKO\tKO definition\n")
                for orf, kos in orf_to_kos.items():
                    for ko in kos:
                        fh.write(f"{orf}\t{ko}\t{ko_to_def.get(ko, '')}\n")

# ================================================================== #
#  ORIGINAL: KOANNOTATION — DISABLED
#  Replaced by kofamscan_direct which runs directly on contigs.faa.
# ================================================================== #

# if SMK_CONFIG:
#     THREADS_KO = 1
# else:
#     THREADS_KO = 32
#
# rule koannotation:
#     input:  expand("{{path}}/temp_splits/Batch_{nb}.faa", nb=range(PRODIGAL_SPLIT))
#     output: out = "{path}/contigs_KEGG_best_hits.tsv"
#     params:
#         root    = "{path}",
#         profile = f"--profile {SMK_CONFIG}" * (SMK_CONFIG != ""),
#     threads: THREADS_KO
#     resources:
#         slurm_partition = get_resource("partition"),
#         mem_mb          = get_resource("mem"),
#     shell: """
#         snakemake -s {LONGFLOW_DIR}/snakenest/ko_annotation.snake \
#             {params.profile} -k --rerun-incomplete \
#             --directory {params.root} --cores {threads} --nolock \
#             --config scripts={SCRIPTS} ROOT={params.root} \
#                 KO_HMM={KO_HMM} KO_HMM_CUTOFFS={KO_HMM_CUTOFFS} \
#                 CONFIG_PATH={params.root}
#         rm -rf {params.root}/.snakemake
#     """

# ================================================================== #
#  ORIGINAL: CAT TAXONOMY  (unchanged — still uses prodigal batches)
# ================================================================== #

if CAT_DB:
    rule CAT_annotation:
        input:
            contigs = "{path}/annotation/temp_splits/Batch_{nb}",
            faa     = "{path}/annotation/temp_splits/Batch_{nb}.faa",
            db      = glob.glob(f"{CAT_DB}/db/*.dmnd"),
        output: ORF2LCA = "{path}/annotation/temp_splits/Batch_{nb}_contigs.contig2classification.txt"
        params: Dir = "{path}/annotation/temp_splits/Batch_{nb}_contigs"
        threads: 32
        resources:
            slurm_partition = get_resource("partition", mult=3, min_size=150000),
            mem_mb          = get_resource("mem", mult=3, min_size=150000),
        shell: """
        if [ -s {input.contigs} ]; then
            {CAT_PATH}/CAT_pack contigs \
                -c {input.contigs} -d {CAT_DB}/db -t {CAT_DB}/tax \
                -p {input.faa} -n {threads} --out_prefix {params.Dir} \
                --top 11 --I_know_what_Im_doing --force \
                --path_to_diamond /hpc-home/kar23heg/bin/diamond
        else
            touch {output}
        fi
        """

    rule CAT_ORF_annotation:
        input:  expand("{{path}}/annotation/temp_splits/Batch_{nb}_contigs.contig2classification.txt", nb=range(PRODIGAL_SPLIT))
        output: "{path}/annotation/CAT_contigs_taxonomy.tsv"
        resources:
            slurm_partition = get_resource("partition"),
            mem_mb          = get_resource("mem"),
        shell: "cat {input} > {output}"

# ================================================================== #
#  ORIGINAL: 16S / BLCA  (unchanged)
# ================================================================== #

rule identify_rRNA:
    input:  "{path}/contigs/contigs.fa"
    output: "{path}/annotation/barrnap_rrna.gff3"
    threads: 32
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb          = get_resource("mem"),
    singularity: SIF_BARRNAP
    shell: "barrnap --threads {threads} {input} > {output}"

localrules: extract_rna_seq

rule extract_rna_seq:
    input:
        contigs = "{path}/contigs/contigs.fa",
        gff     = "{path}/annotation/barrnap_rrna.gff3",
    output:
        seq = "{path}/annotation/16S_seqs.fa",
        bed = "{path}/annotation/16S_bed.tsv",
    run:
        bed = []
        contig_regions = defaultdict(list)
        for line in open(input["gff"]):
            if "16S" in line:
                sline = line.rstrip().split("\t")
                bed.append([sline[0], sline[3], sline[4]])
                contig_regions[sline[0]].append([int(sline[3]), int(sline[4])])
        with open(output["bed"], "w") as handle:
            handle.writelines("%s\n" % "\t".join(line) for line in bed)
        with open(input["contigs"]) as handle, open(output["seq"], "w") as handle_w:
            for header, seq in sfp(handle):
                contig = header.split(" ")[0]
                if contig in contig_regions:
                    for start, end in contig_regions[contig]:
                        handle_w.write(">%s_%s...%s\n%s\n" % (contig, start, end, seq[start:end]))

rule BLCA:
    input:  contigs = "{path}/annotation/16S_seqs.fa"
    output: "{path}/annotation/16S_blca.tsv"
    threads: 32
    params:
        db   = BLCA["db"],
        taxa = BLCA["taxa"],
        path = BLCA["bin"],
    resources:
        slurm_partition = get_resource("partition", min_size=150000),
        mem_mb          = get_resource("mem", min_size=150000),
    shell: """
        module load clustalo/1.2.4
        module load blast+/2.16.0
        cd {params.path}
        python {params.path}/2.blca_main.py \
            -i {input.contigs} -p {threads} -o {output} \
            -r {params.taxa} -q {params.db}
    """

# ================================================================== #
#  STEP 3 — fetchMGs  (parallel SCG route, mOTU 40 universal MGs)
# ================================================================== #

if SIF_FETCHMGS:
    rule fetchMGs:
        input:
            faa = "{path}/annotation/contigs.faa",
            fna = "{path}/annotation/contigs.fna",
        output: done = "{path}/annotation/fetchMGs.done"
        priority: 80
        params: outdir = "{path}/annotation/fetchMGs_out"
        threads: 25
        resources:
            slurm_partition = get_resource("partition"),
            mem_mb          = get_resource("mem"),
        singularity: SIF_FETCHMGS
        shell: r"""
            rm -rf {params.outdir} && mkdir -p {params.outdir}
            fetchMGs extraction {input.faa} gene {params.outdir} \
                -d {input.fna} -t {threads}
            touch {output.done}
        """

# ================================================================== #
#  STEP 4 — GENOMAD  (unchanged)
# ================================================================== #

if GENOMAD_DB:
    rule genomad:
        input:
            contigs = "{path}/contigs/contigs.fa",
            DB      = GENOMAD_DB,
        output: "{path}/annotation/genomad/genomad.done"
        params: output = "{path}/annotation/genomad"
        threads: 25
        resources:
            slurm_partition = get_resource("partition", min_size=150000),
            mem_mb          = get_resource("mem", min_size=150000),
        singularity: SIF_GENOMAD
        shell: """
            genomad end-to-end --cleanup --splits {threads} \
                {input.contigs} {params.output} {GENOMAD_DB} \
                && touch {output}
        """
