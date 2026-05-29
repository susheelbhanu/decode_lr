from os.path import basename, dirname
import glob

# -----------------------------------------------------------------------------
# phyloFlash standalone rules
# -----------------------------------------------------------------------------

DATA = "/ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/metag/results"

PHYLOFLASH_DBHOME = "/ei/.project-scratch/0/0e51ef86-0156-4e79-ad12-c5411c0a5496/databases/phyloflash/138.1"
PHYLOFLASH_SIF = "/ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/decode_lr/singularity_cache/phyloflash_3.3b1.sif"

R1 = {
    basename(dirname(file)): file
    for file in glob.glob(f"{DATA}/preprocessed/reads/*/*_filtered.R1.fq")
}
R1.update({
    basename(dirname(file)): file
    for file in glob.glob(f"{DATA}/preprocessed/reads/*/*_filtered.R1.fq.gz")
})

R2 = {
    basename(dirname(file)): file
    for file in glob.glob(f"{DATA}/preprocessed/reads/*/*_filtered.R2.fq")
}
R2.update({
    basename(dirname(file)): file
    for file in glob.glob(f"{DATA}/preprocessed/reads/*/*_filtered.R2.fq.gz")
})

EXCLUDE_SAMPLES = {"d282"}
PHYLOFLASH_SAMPLES = sorted((set(R1) & set(R2)) - EXCLUDE_SAMPLES)

PHYLOFLASH_DONE = f"{RESULTS_DIR}/phyloflash/{{sample}}/.phyloflash.done"
PHYLOFLASH_TARGETS = expand(PHYLOFLASH_DONE, sample=PHYLOFLASH_SAMPLES)

PHYLOFLASH_HTML = f"{RESULTS_DIR}/phyloflash/{{sample}}/{{sample}}.phyloFlash.html"
PHYLOFLASH_HTML_ABUND = f"{RESULTS_DIR}/phyloflash/{{sample}}/{{sample}}.phyloFlash.html_abundance.tsv"

PHYLOFLASH_HTML_ABUND_TARGETS = expand(PHYLOFLASH_HTML_ABUND, sample=PHYLOFLASH_SAMPLES)

PHYLOFLASH_MERGED_DIR = f"{RESULTS_DIR}/phyloflash/merged"
PHYLOFLASH_MERGED_LONG = f"{PHYLOFLASH_MERGED_DIR}/phyloflash_html_abundance.long.tsv"
PHYLOFLASH_MERGED_WIDE = f"{PHYLOFLASH_MERGED_DIR}/phyloflash_html_abundance.wide.tsv"


rule run_phyloflash:
    input:
        r1=lambda wc: R1[wc.sample],
        r2=lambda wc: R2[wc.sample]
    output:
        done=PHYLOFLASH_DONE
    threads: 24
    resources:
        slurm_partition=get_resource("partition"),
        mem_mb=get_resource("mem")
    params:
        dbhome=PHYLOFLASH_DBHOME,
        sif=PHYLOFLASH_SIF,
        lib=lambda wc: wc.sample,
        outdir=lambda wc: f"{RESULTS_DIR}/phyloflash/{wc.sample}"
    shell:
        r"""
        mkdir -p "{params.outdir}"

        export LC_ALL=C
        export LANG=C

        cd "{params.outdir}"

        singularity exec "{params.sif}" phyloFlash.pl \
          -lib "{params.lib}" -zip -log -taxlevel 6 \
          -read1 "{input.r1}" \
          -read2 "{input.r2}" \
          -CPUs {threads} \
          -dbhome "{params.dbhome}"

        touch "{output.done}"
        """


rule extract_phyloflash_html_abundance:
    input:
        html=PHYLOFLASH_HTML
    output:
        tsv=PHYLOFLASH_HTML_ABUND
    run:
        import re
        import csv
        import html as ihtml
        from pathlib import Path

        in_html = Path(input.html).read_text(encoding="utf-8", errors="ignore")

        # Taxonomy level from the report text
        m_level = re.search(r"Taxonomy summarized at level\s+(\d+)", in_html, flags=re.S)
        tax_level = int(m_level.group(1)) if m_level else None

        # Find the 'Taxonomic affiliation of SSU rRNA reads in library' section
        # and then grab the main abundance table beneath it
        m_block = re.search(
            r'<h3><a href="#" id="taxa-show".*?</h3>\s*<div id="taxa".*?>.*?'
            r'<table class="slimTable">.*?</table>.*?'
            r'<table>\s*(.*?)\s*</table>',
            in_html,
            flags=re.S
        )
        if not m_block:
            raise ValueError(f"Could not find abundance table in {input.html}")

        table_block = m_block.group(1)

        rows = re.findall(
            r"<tr>\s*<td>(.*?)</td>\s*<td>([0-9]+)</td>\s*</tr>",
            table_block,
            flags=re.S
        )

        if not rows:
            raise ValueError(f"No taxon/read rows found in {input.html}")

        Path(output.tsv).parent.mkdir(parents=True, exist_ok=True)

        with open(output.tsv, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["sample", "tax_level", "taxon", "reads"])
            for taxon, reads in rows:
                taxon = ihtml.unescape(re.sub(r"<.*?>", "", taxon)).strip()
                writer.writerow([wildcards.sample, tax_level, taxon, int(reads)])


rule merge_phyloflash_html_abundance:
    input:
        PHYLOFLASH_HTML_ABUND_TARGETS
    output:
        long=PHYLOFLASH_MERGED_LONG,
        wide=PHYLOFLASH_MERGED_WIDE
    run:
        import os
        import pandas as pd

        os.makedirs(os.path.dirname(output.long), exist_ok=True)

        dfs = [pd.read_csv(f, sep="\t") for f in input]
        merged = pd.concat(dfs, ignore_index=True)

        merged = merged.sort_values(["sample", "reads"], ascending=[True, False])
        merged.to_csv(output.long, sep="\t", index=False)

        wide = (
            merged
            .pivot_table(
                index=["tax_level", "taxon"],
                columns="sample",
                values="reads",
                aggfunc="sum",
                fill_value=0
            )
            .reset_index()
        )

        wide.to_csv(output.wide, sep="\t", index=False)


rule all_phyloflash:
    input:
        PHYLOFLASH_TARGETS


rule all_phyloflash_html_abundance:
    input:
        PHYLOFLASH_MERGED_LONG,
        PHYLOFLASH_MERGED_WIDE
