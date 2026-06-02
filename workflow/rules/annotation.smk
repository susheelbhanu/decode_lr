##################################################

# ANNOTATION: LongFlow annotation workflow for
# pre-computed decode_lr assemblies
#
# Adapted from LongFlow (snakenest/annotation.smk)
# Self-contained SLURM resource setup — does NOT include init.smk
# to avoid the full decode_lr sample-loading machinery.
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
                      SLURM_PARTITIONS="", mode="", mult=5, min_size=80000):
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

# SCG detection: use cog_db if present, else checkm HMM
if "cog_db" in ANNOTATION_CFG:
    COG_DB  = ANNOTATION_CFG["cog_db"]
    SCG     = {
        line.rstrip().split("\t")[0]
        for line in open(f"{SCG_DATA}/scg_cogs_min0.97_max1.03_unique_genera.txt")
    }
else:
    CHECKM_DB = ANNOTATION_CFG["checkm"]   # path to checkm database directory
    SCG_HMM   = f"{CHECKM_DB}/hmms/checkm.hmm"
    SCG       = {
        line.rstrip().split("\t")[0]
        for line in open(f"{SCG_DATA}/scg_hmm_selected.txt")
        if line.rstrip().split("\t")[2] == "fine"
    }

NB_SCG = float(len(SCG))

# ------------------------------------------------------------------ #
#  Assembly paths  (flat key -> gz path, defined in config)
# ------------------------------------------------------------------ #
# Keys must not contain '/'.  Example in config:
#   assembly_paths:
#     hifi_individual_b10_d3_0200_con: /ei/.../hifi_individual/b10_d3_0200_con/contigs.fasta.gz
#     hifi_coassembly_bulk: /ei/.../hifi_coassembly_bulk/contigs.fasta.gz
ASSEMBLY_GZ  = config["assembly_paths"]
ASSEMBLIES   = list(ASSEMBLY_GZ.keys())
# print(f"[annotation.smk] ANNOT_DIR:          {ANNOT_DIR}")
# print(f"[annotation.smk] ASSEMBLIES ({len(ASSEMBLIES)}):    {ASSEMBLIES[:3]} ...")
# print(f"[annotation.smk] ANNOTATION_OUTPUT:  (built below)")

# ------------------------------------------------------------------ #
#  Build annotation output list (mirrors LongFlow common.smk logic)
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
# print(f"[annotation.smk] ANNOTATION_OUTPUT:  {ANNOTATION_OUTPUT}")
# print(f"[annotation.smk] Total targets:      {len(ASSEMBLIES) * len(ANNOTATION_OUTPUT)}")

# ------------------------------------------------------------------ #
#  Wildcard constraints
# ------------------------------------------------------------------ #
wildcard_constraints:
    assembly   = "|".join(re.escape(a) for a in ASSEMBLIES),
    annotation = "|".join(DIAMOND.keys()) if DIAMOND else "NO_ANNOTATION",

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
#  PRODIGAL  (parallelised via fasta splits)
# ================================================================== #

rule split_fasta:
    input:  "{annot_dir}/{assembly}/contigs/contigs.fa"
    output: expand("{{annot_dir}}/{{assembly}}/annotation/temp_splits/Batch_{nb}", nb=range(PRODIGAL_SPLIT))
    params: tmp = "{annot_dir}/{assembly}/annotation/temp_splits"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb          = get_resource("mem"),
    singularity: "/hpc-home/kar23heg/singularity_cache/pythonenv_3.9.sif"
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
    singularity: "/hpc-home/kar23heg/singularity_cache/prodigal-gv_2.11.0.sif"
    shell: """
    if [ -s {input} ]; then
        prodigal-gv -i {input} -a {output.faa} -d {output.fna} -f gff -p meta -o {output.gff}
    else
        touch {output}
    fi
    """

rule cat_orfs:
    input:
        faa = expand("{{path}}/temp_splits/Batch_{nb}.faa", nb=range(PRODIGAL_SPLIT)),
        fna = expand("{{path}}/temp_splits/Batch_{nb}.fna", nb=range(PRODIGAL_SPLIT)),
        gff = expand("{{path}}/temp_splits/Batch_{nb}.gff", nb=range(PRODIGAL_SPLIT)),
    output:
        faa = "{path}/contigs.faa",
        fna = "{path}/contigs.fna",
        gff = "{path}/contigs.gff",
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb          = get_resource("mem"),
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
#  SCG ANNOTATION  (COG-db path or CheckM HMM path)
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
        singularity: "/hpc-home/kar23heg/singularity_cache/blast_2.16.0.sif"
        shell: """
            rpsblast -outfmt '6 qseqid sseqid evalue pident length slen qlen' \
                -evalue 0.00001 -query {input} -db {params.db} -out {output} &>{log}
        """

    rule parse_cogs_annotation:
        input:  Batch = expand("{{path}}/annotation/temp_splits/Batch_{nb}.cogs.tsv", nb=range(PRODIGAL_SPLIT))
        output:
            cog = "{path}/annotation/contigs_cogs_best_hits.tsv",
            cat = temp("{path}/annotation/contigs_Cog.out"),
        resources:
            slurm_partition = get_resource("partition"),
            mem_mb          = get_resource("mem"),
        singularity: "/hpc-home/kar23heg/singularity_cache/pythonenv_3.9.sif"
        shell: """
            cat {input} > {output.cat}
            {SCRIPTS}/Filter_Cogs.py {output.cat} \
                --cdd_cog_file {SCG_DATA}/cdd_to_cog.tsv > {output.cog}
        """

else:
    rule hmmsearch:
        input:
            faa = "{path}.faa",
            db  = SCG_HMM,
        output: "{path}_hmm.out"
        log:    "{path}_hmm.log"
        resources:
            slurm_partition = get_resource("partition"),
            mem_mb          = get_resource("mem"),
        singularity: "/hpc-home/kar23heg/singularity_cache/gtdbtk_1.4.0.sif"
        shell: """
        if [ ! -s {input.faa} ]; then
            touch {output}
        else
            hmmsearch --cut_tc --cpu 1 -o /dev/null --noali \
                --domtblout {output} {input.db} {input.faa} >{log} 2>&1
        fi
        """

    rule parse_cogs_annotation:
        input:  Batch = expand("{{path}}/temp_splits/Batch_{nb}_hmm.out", nb=range(PRODIGAL_SPLIT))
        output:
            cat = temp("{path}/contigs_hmm.out"),
            cog = "{path}/contigs_cogs_best_hits.tsv",
        resources:
            slurm_partition = get_resource("partition"),
            mem_mb          = get_resource("mem"),
        singularity: "/hpc-home/kar23heg/singularity_cache/pythonenv_3.9.sif"
        shell: """
            cat {input} > {output.cat}
            {SCRIPTS}/Filter_scg_hmm.py {output.cat} \
                --cog_hmm {SCG_DATA}/scg_hmm_selected.txt {output.cog}
        """

# ================================================================== #
#  SCG EXTRACTION AND CLUSTERING
# ================================================================== #

rule extract_SCG_sequences:
    input:
        annotation = "{filename}_cogs_best_hits.tsv",
        gff        = "{filename}.gff",
        fna        = "{filename}.fna",
    output: "{filename}_SCG.fna"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb          = get_resource("mem"),
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
    singularity: "/hpc-home/kar23heg/singularity_cache/mmseqs2_latest.sif"
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
    singularity: "/hpc-home/kar23heg/singularity_cache/mmseqs2_latest.sif"
    shell: """
        mmseqs easy-cluster {input} {params.out} {params.tmp} \
            --min-seq-id {wildcards.pid} -c 0.8 \
            --cov-mode 1 --alignment-mode 3 --threads {threads}
    """

# ================================================================== #
#  CONTIG QUALITY (bed + summary)
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

        get_cogs     = lambda x: Counter([el.split("_")[0] for el in contigs_to_scgs[x]])
        completion   = lambda x: sum(v >= 1 for v in get_cogs(x).values()) / NB_SCG
        contamination= lambda x: (sum(get_cogs(x).values()) - len(get_cogs(x))) / NB_SCG

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
#  DIAMOND ANNOTATION
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
    singularity: "/hpc-home/kar23heg/singularity_cache/diamond_2.0.9.sif"
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
        annotation     = lambda w: DIAMOND[w.annotation]["annotation"],
        Bitscore       = lambda w: DIAMOND[w.annotation]["filter"][0],
        Evalue         = lambda w: DIAMOND[w.annotation]["filter"][1],
        PID            = lambda w: DIAMOND[w.annotation]["filter"][2],
        subject_pid    = lambda w: DIAMOND[w.annotation]["filter"][3],
        subject_cov    = lambda w: DIAMOND[w.annotation]["filter"][4],
        query_cov      = lambda w: DIAMOND[w.annotation]["filter"][5],
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
#  KO ANNOTATION  (kofamscan, nested snakemake)
# ================================================================== #

if SMK_CONFIG:
    THREADS_KO = 1
else:
    THREADS_KO = 32

rule koannotation:
    input:  expand("{{path}}/temp_splits/Batch_{nb}.faa", nb=range(PRODIGAL_SPLIT))
    output: out = "{path}/contigs_KEGG_best_hits.tsv"
    params:
        root    = "{path}",
        profile = f"--profile {SMK_CONFIG}" * (SMK_CONFIG != ""),
    threads: THREADS_KO
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb          = get_resource("mem"),
    shell: """
        snakemake -s {LONGFLOW_DIR}/snakenest/ko_annotation.snake \
            {params.profile} -k --rerun-incomplete \
            --directory {params.root} --cores {threads} --nolock \
            --config scripts={SCRIPTS} ROOT={params.root} \
                KO_HMM={KO_HMM} KO_HMM_CUTOFFS={KO_HMM_CUTOFFS} \
                CONFIG_PATH={params.root}
        rm -rf {params.root}/.snakemake
    """

# ================================================================== #
#  CAT TAXONOMY
# ================================================================== #

if CAT_DB:
    rule CAT_annotation:
        input:
            contigs = "{path}/annotation/temp_splits/Batch_{nb}",
            faa     = "{path}/annotation/temp_splits/Batch_{nb}.faa",
            db      = glob.glob(f"{CAT_DB}/db/*.dmnd"),
        output: ORF2LCA = "{path}/annotation/temp_splits/Batch_{nb}_contigs.contig2classification.txt"
        params: Dir = "{path}/annotation/temp_splits/Batch_{nb}_contigs"
        threads: 48
        resources:
            slurm_partition = get_resource("partition", mult=7),
            mem_mb          = get_resource("mem", mult=7),
        shell: """
        if [ -s {input.contigs} ]; then
#            module load diamond/2.1.15.169
            {CAT_PATH}/CAT_pack contigs \
                -c {input.contigs} -d {CAT_DB}/db -t {CAT_DB}/tax \
                -p {input.faa} -n {threads} --out_prefix {params.Dir} \
                --top 11 --I_know_what_Im_doing --force --path_to_diamond /hpc-home/kar23heg/bin/diamond
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
#  16S / BLCA
# ================================================================== #

rule identify_rRNA:
    input:  "{path}/contigs/contigs.fa"
    output: "{path}/annotation/barrnap_rrna.gff3"
    threads: 32
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb          = get_resource("mem"),
    singularity: "/hpc-home/kar23heg/singularity_cache/barrnap_0.9.sif"
    shell: "barrnap --threads {threads} {input} > {output}"

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
        slurm_partition = get_resource("partition"),
        mem_mb          = get_resource("mem"),
    shell: """
        module load clustalo/1.2.4 && module load blast+/2.16.0 && cd {params.path}
        python {params.path}/2.blca_main.py \
            -i {input.contigs} -p {threads} -o {output} \
            -r {params.taxa} -q {params.db}
    """

# ================================================================== #
#  GENOMAD
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
            slurm_partition = get_resource("partition"),
            mem_mb          = get_resource("mem"),
        singularity: "/hpc-home/kar23heg/singularity_cache/genomad_1.8.0.sif"
        shell: """
            genomad end-to-end --cleanup --splits {threads} \
                {input.contigs} {params.output} {GENOMAD_DB} \
                && touch {output}
        """

