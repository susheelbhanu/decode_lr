#############################################
# workflow/Snakefile_binning_standalone.smk
#############################################

import os
import glob
import yaml
from pathlib import Path

with open(config["lr_config"]) as fh:
    lr_cfg = yaml.safe_load(fh)
RESULTS_DIR = lr_cfg["results_dir"]

RUN = config.get("run_name", "binning_run")
OUTDIR = f"{RESULTS_DIR}/binning/{RUN}"

ASSEMBLY_GZ = config["assembly"]
LR_READS = config["lr_reads"]  # dict sample -> fastq.gz

USE_SR = bool(config.get("use_sr", False))
SR_LIST = config.get("sr_sample_list")
SR_ROOT = config.get("sr_data_root")
SR_R1_GLOB = config.get("sr_r1_glob", "{root}/preprocessed/reads/*/*_filtered.R1.fq")
SR_R2_GLOB = config.get("sr_r2_glob", "{root}/preprocessed/reads/*/*_filtered.R2.fq")

BINNERS = config.get("binners", ["metabat2", "semibin2", "concoct"])
MM2_PRESET = config.get("minimap2_preset", "-ax map-hifi")
THREADS_MAP = int(config.get("threads_map", 32))
THREADS_BIN = int(config.get("threads_bin", 24))
TMPDIR = config.get("tmpdir", "/tmp")

MB2 = config.get("metabat2", {})
MB2_MINLEN = int(MB2.get("min_contig_len", 2000))
MB2_OPTS = str(MB2.get("extra_opts", "")).strip()

CC = config.get("concoct", {})
CC_MINLEN = int(CC.get("min_contig_len", 1000))
CC_CHUNK = int(CC.get("chunk_len", 10000))
CC_OVERLAP = int(CC.get("overlap_len", 0))

CHECKM2_THREADS = int(config.get("checkm2", {}).get("threads", THREADS_BIN))
GTDB_THREADS = int(config.get("gtdbtk", {}).get("threads", THREADS_BIN))

Path(OUTDIR).mkdir(parents=True, exist_ok=True)
Path(TMPDIR).mkdir(parents=True, exist_ok=True)

def load_sr_samples():
    if not USE_SR:
        return []
    if not SR_LIST:
        raise ValueError("use_sr=true but sr_sample_list not set")
    with open(SR_LIST) as fh:
        return [x.strip() for x in fh if x.strip()]

def discover_sr_pairs():
    if not USE_SR:
        return {}
    if not SR_ROOT:
        raise ValueError("use_sr=true but sr_data_root not set")

    r1_files = glob.glob(SR_R1_GLOB.format(root=SR_ROOT))
    r2_files = glob.glob(SR_R2_GLOB.format(root=SR_ROOT))

    r1 = {os.path.basename(os.path.dirname(f)): f for f in r1_files}
    r2 = {os.path.basename(os.path.dirname(f)): f for f in r2_files}

    wanted = set(load_sr_samples())
    pairs = {}
    missing = []
    for s in wanted:
        if s in r1 and s in r2:
            pairs[s] = (r1[s], r2[s])
        else:
            missing.append(s)

    if missing:
        print(f"[WARN] {len(missing)} SR samples listed but not found via glob (first 10): {missing[:10]}")
    return pairs

SR_PAIRS = discover_sr_pairs() if USE_SR else {}
SR_SAMPLES = sorted(SR_PAIRS.keys())

# ---------- paths ----------
ASM_FA = f"{OUTDIR}/assembly/contigs.fa"
ASM_FAI = f"{ASM_FA}.fai"

LR_BAMS = expand(f"{OUTDIR}/mapping/LR/{{sample}}.sorted.bam", sample=list(LR_READS.keys()))
LR_MERGED = f"{OUTDIR}/mapping/LR/all.sorted.bam"
LR_MERGED_BAI = f"{LR_MERGED}.bai"

SR_BAMS = expand(f"{OUTDIR}/mapping/SR/{{sample}}.sorted.bam", sample=SR_SAMPLES) if USE_SR else []
SR_MERGED = f"{OUTDIR}/mapping/SR/rhizo.sorted.bam"
SR_MERGED_BAI = f"{SR_MERGED}.bai"

DEPTH = f"{OUTDIR}/depth/depth.txt"

MB2_DIR = f"{OUTDIR}/binning/metabat2/bins"
MB2_DONE = f"{OUTDIR}/binning/metabat2/.metabat2.done"

SB2_DIR = f"{OUTDIR}/binning/semibin2"
SB2_DONE = f"{OUTDIR}/binning/semibin2/.semibin2.done"

CC_DIR = f"{OUTDIR}/concoct"
CC_FILTER_FA = f"{CC_DIR}/contigs.min{CC_MINLEN}.fa"
CC_CUT_FA = f"{CC_DIR}/contigs.cut.fa"
CC_BED = f"{CC_DIR}/contigs.cut.bed"
CC_COV = f"{CC_DIR}/coverage.tsv"
CC_CLUSTER = f"{CC_DIR}/clustering.csv"
CC_BINS_DIR = f"{OUTDIR}/binning/concoct/bins"
CC_DONE = f"{OUTDIR}/binning/concoct/.concoct.done"

CHECKM2_OUT = f"{OUTDIR}/qc/checkm2"
CHECKM2_TSV = f"{CHECKM2_OUT}/quality_report.tsv"

GTDB_OUT = f"{OUTDIR}/qc/gtdbtk"
GTDB_SUM = f"{GTDB_OUT}/gtdbtk.bac120.summary.tsv"

TREE = f"{OUTDIR}/qc/tree/fasttree.nwk"
BINS_LIST = f"{OUTDIR}/qc/bins.list"

ALL_DONE = f"{OUTDIR}/.binning.done"

# ---------- containers (matching your binning.smk where possible) ----------
CACHE = "/ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/decode_lr/singularity_cache"

IMG_BWA_SAMTOOLS = f"{CACHE}/bwasamtools_1.10.sif"
IMG_SEMIBIN      = f"{CACHE}/semibin2.sif"
IMG_CHECKM2      = f"{CACHE}/checkm2.sif"
IMG_GTDBTK       = f"{CACHE}/gtdbtk_2.4.0.sif"
IMG_FASTTREE     = f"{CACHE}/fasttree_2.1.11.sif"

# Needed extras for standalone mapping/binning utilities
IMG_MINIMAP2     = f"{CACHE}/minimap2_2.28.sif"
IMG_PIGZ         = f"{CACHE}/pigz_2.4.sif"
IMG_SAMTOOLS     = f"{CACHE}/samtools_1.19.2.sif"
IMG_METABAT2     = f"{CACHE}/metabat2_2.15.sif"
IMG_CONCOCT      = f"{CACHE}/concoct_1.1.0.sif"
IMG_SEQKIT       = f"{CACHE}/seqkit_2.6.1.sif"
IMG_MINIMAP2_SAMTOOLS = "docker://quay.io/biocontainers/mulled-v2-66534bcbb703196e9da7354972e259e875080e7e:afc8599426f4b66471e4eb441235123d47228d71-0"

# ---------- targets ----------
targets = []
if "metabat2" in BINNERS:
    targets.append(MB2_DONE)
if "semibin2" in BINNERS:
    targets.append(SB2_DONE)
if "concoct" in BINNERS:
    targets.append(CC_DONE)
targets += [CHECKM2_TSV, GTDB_SUM, TREE, ALL_DONE]

rule all:
    input: targets

# ---------- assembly ----------
rule prepare_assembly:
    input:
        fa_gz=lambda wc: config["assembly"]
    output:
        fa=ASM_FA
    threads: 1
    shell:
        r"""
        mkdir -p "$(dirname "{output.fa}")"
        if [[ "{input.fa_gz}" == *.gz ]]; then
          gzip -dc "{input.fa_gz}" > "{output.fa}"
        else
          ln -sf "{input.fa_gz}" "{output.fa}"
        fi
        """

rule faidx_assembly:
    input: ASM_FA
    output: ASM_FAI
    singularity: IMG_SAMTOOLS
    shell: r"samtools faidx {input}"

# ---------- LR mapping per sample ----------
rule map_lr_sample:
    input:
        fa=ASM_FA,
        fai=ASM_FAI,
        fq=lambda wc: LR_READS[wc.sample]
    output:
        bam=f"{OUTDIR}/mapping/LR/{{sample}}.sorted.bam",
        bai=f"{OUTDIR}/mapping/LR/{{sample}}.sorted.bam.bai"
    threads: THREADS_MAP
    singularity: IMG_MINIMAP2_SAMTOOLS # IMG_MINIMAP2
    shell:
        r"""
        mkdir -p {OUTDIR}/mapping/LR
        minimap2 {MM2_PRESET} -t {threads} "{input.fa}" "{input.fq}" \
          | samtools view -@ {threads} -b - \
          | samtools sort -@ {threads} -T "{TMPDIR}/mm2.{wildcards.sample}" -o "{output.bam}"
        samtools index "{output.bam}"
        """

rule merge_lr_bams:
    input: LR_BAMS
    output:
        bam=LR_MERGED,
        bai=LR_MERGED_BAI
    threads: 8
    singularity: IMG_SAMTOOLS
    shell:
        r"""
        mkdir -p {OUTDIR}/mapping/LR
        samtools merge -@ {threads} -f "{output.bam}" {input}
        samtools index "{output.bam}"
        """

# ---------- SR mapping per sample (optional) ----------
rule map_sr_sample:
    input:
        fa=ASM_FA,
        fai=ASM_FAI,
        r1=lambda wc: SR_PAIRS[wc.sample][0],
        r2=lambda wc: SR_PAIRS[wc.sample][1]
    output:
        bam=f"{OUTDIR}/mapping/SR/{{sample}}.sorted.bam",
        bai=f"{OUTDIR}/mapping/SR/{{sample}}.sorted.bam.bai"
    threads: THREADS_MAP
    singularity: IMG_BWA_SAMTOOLS
    shell:
        r"""
        mkdir -p {OUTDIR}/mapping/SR
        bwa index "{input.fa}" 1>/dev/null 2>/dev/null || true
        bwa mem -t {threads} "{input.fa}" "{input.r1}" "{input.r2}" \
          | samtools view -@ {threads} -b - \
          | samtools sort -@ {threads} -T "{TMPDIR}/bwa.{wildcards.sample}" -o "{output.bam}"
        samtools index "{output.bam}"
        """

rule merge_sr_bams:
    input: SR_BAMS
    output:
        bam=SR_MERGED,
        bai=SR_MERGED_BAI
    threads: 8
    singularity: IMG_SAMTOOLS
    shell:
        r"""
        mkdir -p {OUTDIR}/mapping/SR

        # If there are no SR BAMs (e.g. none discovered), create empty placeholders
        if [[ -z "{input}" ]]; then
          touch "{output.bam}" "{output.bai}"
          exit 0
        fi

        samtools merge -@ {threads} -f "{output.bam}" {input}
        samtools index "{output.bam}"
        """

# ---------- depth ----------
rule depth:
    input:
        lr=LR_MERGED,
        sr=SR_MERGED
    output: DEPTH
    singularity: IMG_METABAT2
    shell:
        r"""
        mkdir -p {OUTDIR}/depth
        if [[ -s "{input.sr}" ]]; then
          jgi_summarize_bam_contig_depths --minContigLength {MB2_MINLEN} --outputDepth "{output}" "{input.lr}" "{input.sr}"
        else
          jgi_summarize_bam_contig_depths --minContigLength {MB2_MINLEN} --outputDepth "{output}" "{input.lr}"
        fi
        """

# ---------- MetaBAT2 ----------
rule metabat2:
    input:
        fa=ASM_FA,
        depth=DEPTH
    output: MB2_DONE
    threads: THREADS_BIN
    singularity: IMG_METABAT2
    shell:
        r"""
        mkdir -p "{MB2_DIR}"
        metabat2 -i "{input.fa}" -a "{input.depth}" -o "{MB2_DIR}/bin" -t {threads} {MB2_OPTS}

        i=1
        for f in {MB2_DIR}/bin.*.fa; do
          [[ -e "$f" ]] || continue
          printf -v nn "%03d" "$i"
          mv "$f" "{MB2_DIR}/Bin_${{nn}}.fa"
          i=$((i+1))
        done

        touch "{output}"
        """

# ---------- SemiBin2 ----------
rule semibin2:
    input:
        fa=ASM_FA,
        bam=LR_MERGED
    output: SB2_DONE
    threads: THREADS_MAP
    singularity: IMG_SEMIBIN
    shell:
        r"""
        mkdir -p "{SB2_DIR}"
        SemiBin2 single_easy_bin \
          --contig_fasta "{input.fa}" \
          --bam "{input.bam}" \
          --output "{SB2_DIR}" \
          --threads {threads}
        touch "{output}"
        """

# ---------- CONCOCT ----------
rule concoct_filter_contigs:
    input: ASM_FA
    output: CC_FILTER_FA
    singularity: IMG_SEQKIT
    shell:
        r"""
        mkdir -p "{CC_DIR}"
        seqkit seq -m {CC_MINLEN} "{input}" > "{output}"
        """

rule concoct_cutup:
    input: CC_FILTER_FA
    output:
        cut=CC_CUT_FA,
        bed=CC_BED
    singularity: IMG_CONCOCT
    shell:
        r"""
        mkdir -p "{CC_DIR}"
        cut_up_fasta.py "{input}" -c {CC_CHUNK} -o {CC_OVERLAP} --merge_last > "{output.cut}"
        cut_up_fasta.py "{input}" -c {CC_CHUNK} -o {CC_OVERLAP} --merge_last --bedfile > "{output.bed}"
        """

rule concoct_coverage:
    input:
        bed=CC_BED,
        lr=LR_MERGED,
        sr=SR_MERGED
    output: CC_COV
    singularity: IMG_CONCOCT
    shell:
        r"""
        mkdir -p "{CC_DIR}"
        if [[ -s "{input.sr}" ]]; then
          concoct_coverage_table.py "{input.bed}" "{input.lr}" "{input.sr}" > "{output}"
        else
          concoct_coverage_table.py "{input.bed}" "{input.lr}" > "{output}"
        fi
        """

rule concoct_cluster:
    input:
        cut=CC_CUT_FA,
        cov=CC_COV
    output: CC_CLUSTER
    threads: THREADS_BIN
    singularity: IMG_CONCOCT
    shell:
        r"""
        mkdir -p "{CC_DIR}"
        concoct --composition_file "{input.cut}" --coverage_file "{input.cov}" -b "{CC_DIR}/" -t {threads}
        cp "{CC_DIR}/clustering.csv" "{output}"
        """

rule concoct_extract_bins:
    input:
        fa=CC_FILTER_FA,
        clustering=CC_CLUSTER
    output: CC_DONE
    singularity: IMG_CONCOCT
    shell:
        r"""
        mkdir -p "{CC_BINS_DIR}"
        extract_fasta_bins.py "{input.fa}" "{input.clustering}" --output_path "{CC_BINS_DIR}"
        touch "{output}"
        """

# ---------- bins list for QC ----------
rule make_bins_list:
    input:
        mb2=MB2_DONE if "metabat2" in BINNERS else [],
        sb2=SB2_DONE if "semibin2" in BINNERS else [],
        cc=CC_DONE if "concoct" in BINNERS else []
    output: BINS_LIST
    shell:
        r"""
        mkdir -p "{OUTDIR}/qc"
        python - << 'PY'
import glob
bins=[]
bins += glob.glob(r"{MB2_DIR}/Bin_*.fa")
bins += glob.glob(r"{SB2_DIR}/output_bins/*.fa")
bins += glob.glob(r"{CC_BINS_DIR}/*.fa")
bins = sorted(set(bins))
with open(r"{BINS_LIST}","w") as fh:
    for b in bins:
        fh.write(b+"\n")
print("bins:", len(bins))
PY
        """

# ---------- CheckM2 ----------
rule checkm2:
    input: BINS_LIST
    output: CHECKM2_TSV
    threads: CHECKM2_THREADS
    singularity: IMG_CHECKM2
    shell:
        r"""
        mkdir -p "{CHECKM2_OUT}"
        binsdir="{CHECKM2_OUT}/bins"
        mkdir -p "$binsdir"
        while read -r f; do
          ln -sf "$f" "$binsdir/$(basename "$f")"
        done < "{input}"
        checkm2 predict --threads {threads} --input "$binsdir" --output-directory "{CHECKM2_OUT}"
        test -s "{output}"
        """

# ---------- GTDB-Tk ----------
rule gtdbtk:
    input: BINS_LIST
    output: GTDB_SUM
    threads: GTDB_THREADS
    singularity: IMG_GTDBTK
    shell:
        r"""
        mkdir -p "{GTDB_OUT}"
        binsdir="{GTDB_OUT}/bins"
        mkdir -p "$binsdir"
        while read -r f; do
          ln -sf "$f" "$binsdir/$(basename "$f")"
        done < "{input}"
        gtdbtk classify_wf --genome_dir "$binsdir" --out_dir "{GTDB_OUT}" --cpus {threads} --skip_ani_screen
        test -s "{output}"
        """

# ---------- FastTree ----------
rule fasttree:
    input: GTDB_SUM
    output: TREE
    singularity: IMG_FASTTREE
    shell:
        r"""
        mkdir -p "{OUTDIR}/qc/tree"
        msa="{GTDB_OUT}/align/gtdbtk.bac120.msa.fasta"
        if [[ -s "$msa" ]]; then
          FastTree -nt "$msa" > "{output}"
        else
          echo "();" > "{output}"
        fi
        """

# ---------- Final marker ----------
rule binning_done:
    input:
        mb2=MB2_DONE if "metabat2" in BINNERS else [],
        sb2=SB2_DONE if "semibin2" in BINNERS else [],
        cc=CC_DONE if "concoct" in BINNERS else [],
        checkm2=CHECKM2_TSV,
        gtdb=GTDB_SUM,
        tree=TREE
    output: ALL_DONE
    shell: r'touch "{output}"'
