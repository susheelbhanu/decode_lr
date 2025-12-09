import sys
sys.path.append('/ei/.project-scratch/0/0e51ef86-0156-4e79-ad12-c5411c0a5496/seb/repos/Snakescripts/python')

import os 
import glob
import numpy as np
from os.path import basename,dirname,realpath, isfile
from collections import defaultdict, Counter
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
from utils import get_resource, get_resource_real, matrix_write, load_matrix, get_mag_dmag

SCRIPTS = "/ei/.project-scratch/0/0e51ef86-0156-4e79-ad12-c5411c0a5496/seb/repos/Metahood/scripts"
SCRIPTS2 = "/ei/.project-scratch/0/0e51ef86-0156-4e79-ad12-c5411c0a5496/seb/repos/LongFlow/snakenest/python"
SCG_DATA = "/ei/.project-scratch/0/0e51ef86-0156-4e79-ad12-c5411c0a5496/seb/repos/Metahood/scg_data"
COG_DB = "/ei/.project-scratch/0/0e51ef86-0156-4e79-ad12-c5411c0a5496/databases/rpsblast_cog_db/Cog"

CONTAINER_DIR = "/ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/decode_lr/containers"

BWASAMTOOLS_IMG = f"{CONTAINER_DIR}/bwasamtools_1.10.sif"
BEDTOOLS_IMG    = f"{CONTAINER_DIR}/bedtools_2.29.2.sif"
PRODIGAL_IMG    = f"{CONTAINER_DIR}/prodigal-gv_2.11.0.sif"
BLAST_IMG       = f"{CONTAINER_DIR}/blast_2.16.0.sif"
PYTHONENV_IMG   = f"{CONTAINER_DIR}/pythonenv_3.9.sif"
SYLPH_IMG       = f"{CONTAINER_DIR}/sylph_0.9.0.sif"

ROOT = "/ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7/CEHsoil/HiFi/assemblies/mags_with_rhyzo"
DATA = "/ei/projects/5/542de014-1e71-4955-945a-5d2ab09567a7/CEH_soil_project/metag/results"
DREP95 = "/ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7/CEHsoil/HiFi/assemblies/mags_with_rhyzo/drep"
SUMMARY = "/ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7/CEHsoil/HiFi/assemblies/mags_with_rhyzo/results/contig_MAGs_dMAGs.tsv"
OUT = f'{ROOT}/MAGs/profile/scg_95'

SPIKIN_FOLD = '/ei/.project-scratch/0/0e51ef86-0156-4e79-ad12-c5411c0a5496/databases/spike_in_refs'
SPIKE_IN = {'Allobacillus_halotolerans', 'Imtechella_halotolerans'}


MAGs = {realpath(f) for f in glob.glob(f"{ROOT}/all_mags/*/*.fa")}
MAGs|={f"{SPIKIN_FOLD}/{spk}.fasta" for spk in SPIKE_IN}
DMAGs = {realpath(f) for f in glob.glob(f"{DREP95}/dereplicated_genomes/*.fa")}
DMAGs|= {f"{SPIKIN_FOLD}/{spk}.fasta" for spk in SPIKE_IN}

R1 = {basename(dirname(file)):file for file in glob.glob(f"{DATA}/preprocessed/reads/*/*_filtered.R1.fq")}
# R1.update({basename(dirname(file)):file for file in glob.glob(f"{DATA}/*SPIKE/*R1.fastq.gz")})
R2 = {basename(dirname(file)):file for file in glob.glob(f"{DATA}/preprocessed/reads/*/*_filtered.R2.fq")}
# R2.update({basename(dirname(file)):file for file in glob.glob(f"{DATA}/*SPIKE/*R2.fastq.gz")})

# ASM_folder_SR = {"COA_SR":f"{ROOT}/WGS/assemblies/coassembly"}
# ASM_folder_SR.update({"%s_SR"%basename(file):file for file in glob.glob(f"{ROOT}/WGS/assemblies/per_mice/*")})
# ASM_folder_LR = {"COA_LR":f"/ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7/CEHsoil/HiFi/assemblies/COA/metamdbg"}
#ASM_folder_LR.update({"%s_LR"%basename(dirname(file)):file for file in glob.glob(f"{ROOT}/HiFi/assemblies/SSA/*/metamdbg")})
ASM_folder_LR = {"%s_LR"%basename(dirname(file)):file for file in glob.glob(f"/ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7/CEHsoil/HiFi/assemblies/SSA/*/metamdbg")}

# ASM_folder = ASM_folder_SR
#ASM_folder.update(ASM_folder_LR)
#for ignore in {'EbNeg_SR', 'ExNeg_SR'}:
#    del ASM_folder[ignore]
ASM_folder = ASM_folder_LR

################  Read mapping ################
# get orfs collections for each mag
# for now just get all the contigs of all the scg of each 95% dmag
# build the bed file as well

rule results:
    input: f"{OUT}/drep_95.cov",
           f"{OUT}/drep_95_orfs.cov",
           f"{OUT}/drep_95_orfs_detail.cov",
           [f"{OUT}/../sylph/cov/{sample}.tsv" for sample in R1]

rule create_orf_database:
#    input: expand("%s/{spk}_SCG.fna"%SPIKIN_FOLD,spk=SPIKE_IN)
    output: fa = "{path}/dmag_scg.fa",
            bed = "{path}/dmag_scg.bed"
    run:
        mag_to_dmag,dmag_to_mags = get_mag_dmag(DREP95)
        asm_mags = defaultdict(set)
        for mag in mag_to_dmag:
            if mag not in SPIKE_IN:
                asm,mag = mag.split("_Bin")
                mag = f"Bin_{mag}"
                asm_mags[asm].add(mag)
        sanity_check = set()
        with open(output["fa"],"w") as handle_fa, open(output["bed"],"w") as handle_bed:
            for asm,fold in ASM_folder.items():

                # get scg
                scg_def_file = f"{fold}/annotation/contigs_SCG.fna"
                scg_def = {header.split(" ")[0]:header.split(" ")[1] for header,seq in sfp(open(scg_def_file))}

                # get bed
                orfs_def = f"{fold}/annotation/orf.bed"
                orf_scg = {line.rstrip().split("\t")[3]:line.rstrip().split("\t") for line in open(orfs_def) if line.rstrip().split("\t")[3] in scg_def}

                contigs ={contig for contig,start,end,orf in orf_scg.values()}

                # get contigs
                contig_mag = {}
                if "LR" in asm:
                    clustering = f"{fold}/binning/consensus_LR/clustering_consensus_LR.csv"
                    LR_mags = {basename(f).replace(".fa",""):f for f in glob.glob(f"{fold}/MAGs/mags_LRSR/*.fa")}
                    for mag,mag_path in LR_mags.items():
                        for header,seq in sfp(open(mag_path)):
                            ctg = header.split(" ")[0]
                            if ctg in contigs:
                                mag = f"{mag}"
                                mag_asm = f"{asm}_{mag}"
                                contig_mag[ctg] = mag_asm
                else:
                    clustering = f"{fold}/binning/consensus/clustering_consensus.csv"
                for line in open(clustering):
                    ctg,mag = line.rstrip().split(",")
                    if ctg in contigs:
                        mag = f"Bin_{mag}"
                        mag_asm = f"{asm}_{mag}"
                        contig_mag[ctg] = mag_asm
                contigs = set(contig_mag.keys())
                orf_scg = {key:value for key,value in orf_scg.items() if value[0] in contigs}

                # write to output
                sorted_contig = []
                asm_file = f"{fold}/contigs/contigs.fa"
                for header,seq in sfp(open(asm_file)):
                    header = header.split(" ")[0]
                    if header in contigs:
                        mag = contig_mag[header]
                        handle_fa.write(f">{mag}__{header}\n{seq}\n")
                        sorted_contig.append(f"{mag}__{header}")

                # sort everything for bedtool
                contig_beds = defaultdict(list)
                for orf,info in orf_scg.items():
                    contig = info[0]
                    mag = contig_mag[contig]
                    mcontig = f"{mag}__{contig}"
                    if contig in contigs:
                        info[0] = mcontig
                        orf = info[3]
                        morf = f"{mag}__{orf}"
                        info[3] = morf
                        infos = info+[scg_def[orf],contig_mag[contig]]
                        contig_beds[mcontig].append(infos)

                contig_beds = {ctg:sorted(bed,key=lambda x:int(x[1])) for ctg,bed in contig_beds.items()}
                for contig in sorted_contig:
                    handle_bed.writelines("%s\n"%"\t".join(bed) for bed in contig_beds[contig])
                sanity_check|=set(contig_mag.values())
        print(len(sanity_check))

rule deal_with_spike_in:
    input: spk = expand("%s/{spk}_SCG.fna"%SPIKIN_FOLD,spk=SPIKE_IN),
           spk_bed = expand("%s/{spk}.bed"%SPIKIN_FOLD,spk=SPIKE_IN),
           fa = "{path}/dmag_scg.fa",
           bed = "{path}/dmag_scg.bed"
    output: fa = "{path}/dmag_spk_scg.fa",
            bed = "{path}/dmag_spk_scg.bed"
    run: 
        with open(output["fa"],"w") as handle_fa, open(output["bed"],"w") as handle_bed:
            for spk in SPIKE_IN:

                # scg
                scg_def = {header.split(" ")[0]:header.split(" ")[1] for header,seq in sfp(open(f"{SPIKIN_FOLD}/{spk}_SCG.fna"))}
                contigs = {"_".join(orf.split("_")[:-1]) for orf in scg_def}

                sorted_contig = []
                # write fasta
                for header,seq in sfp(open(f"{SPIKIN_FOLD}/{spk}.fasta")):
                    header = header.split(" ")[0]
                    if header in contigs:
                        handle_fa.write(f">{spk}__{header}\n%s{seq}\n")
                        sorted_contig.append(f"{spk}__{header}")

                # bed
                orfs_def = f"{SPIKIN_FOLD}/{spk}.bed"
                orf_scg = {line.rstrip().split("\t")[3]:line.rstrip().split("\t") for line in open(orfs_def) if line.rstrip().split("\t")[3] in scg_def}
                
                contig_beds = defaultdict(list)
                for orf,info in orf_scg.items():
                    contig = info[0]
                    mcontig = f"{spk}__{contig}"
                    info[0] = mcontig
                    orf = info[3]
                    morf = f"{spk}__{orf}"
                    info[3] = morf
                    infos = info+[scg_def[orf],spk]
                    contig_beds[mcontig].append(infos)
                contig_beds = {ctg:sorted(bed,key=lambda x:int(x[1])) for ctg,bed in contig_beds.items()}

                # write bed
                for contig in sorted_contig:
                    handle_bed.writelines("%s\n"%"\t".join(bed) for bed in contig_beds[contig])

        # concat:
        os.system("cat %s >> %s"%(input["fa"],output["fa"]))
        os.system("cat %s >> %s"%(input["bed"],output["bed"]))




rule bwa_index:
    input:   "{path}/dmag_spk_scg.fa",
    output:  touch("{path}/index.done")
    log:     "{path}/index.log"
    params : 1000000
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb=get_resource("mem"),
    message: "Building bwa index for {input}"
    singularity: BWASAMTOOLS_IMG
    shell:   "bwa index -b {params} {input} &> {log}"

rule bwa_mem_to_bam:
    input:   index="{group}/index.done",
             contigs="{group}/dmag_spk_scg.fa",
             R1 = lambda w:R1[w.sample],
             R2 = lambda w:R2[w.sample]
    output:  "{group}/bam/{sample}_mapped_sorted.bam"
    threads: 4
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")
    log:     "{group}/bam/{sample}_map.log"
    singularity: BWASAMTOOLS_IMG
    shell:   "bwa mem -t {threads} {input.contigs} {input.R1} {input.R2} 2>{log} | samtools view  -b -F 4 -@{threads} - | samtools sort -@{threads} - > {output}"

rule index:
    input: bam = "{path}_mapped_sorted.bam"
    output: "{path}_mapped_sorted.bam.bai"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")
    singularity: BWASAMTOOLS_IMG
    shell: "samtools index {input.bam}"

rule spike_in_def:
    output: "{path}/spike_in_def.tsv"
    run:
        with open(output[0],"w") as handle:
            handle.writelines(f"{spk}\n" for spk in SPIKE_IN)

rule drep_clu_scg:
    input: bam ="{path}/bam/{sample}_mapped_sorted.bam",
           bai ="{path}/bam/{sample}_mapped_sorted.bam.bai",
           bed = "{path}/dmag_spk_scg.bed",
           spk = "{path}/spike_in_def.tsv"
    output: orf = "{path}/cov/{sample}_drep_detail.cov",
            cov = "{path}/cov/{sample}_drep.cov",
            bam_filt = "{path}/bam/{sample}_filtered.bam"
    params: derep = SUMMARY,
            prefix = "{path}/cov/{sample}",
    threads: 4
    resources:
        slurm_partition = get_resource("partition",mult=10),
        mem_mb = get_resource("mem",mult=10)
    shell: "{SCRIPTS2}/drep_clu_scg.py {input.bam} {params.derep} {input.bed} {input.spk} {params.prefix}"


rule concat_results_drep:
    input: cov = expand("{{path}}/cov/{sample}_drep.cov",sample=R1),
           table = expand("{{path}}/cov/{sample}_drep_detail.cov",sample=R1)
    output: cov = '{path}/drep_95.cov',
            table = '{path}/drep_95_detail.cov'
    run:
        sorted_mags = [line.rstrip().split("\t")[0] for index,line in enumerate(open(input[0])) if index>0]
        sorted_samples = sorted(R1)
        matcov = np.zeros((len(sorted_mags),len(sorted_samples)))
        for file in input["cov"]:
            with open(file) as handle:
                sample = next(handle).rstrip().split("\t")[1]
                index_col = sorted_samples.index(sample)
                for index_row,line in enumerate(handle):
                    matcov[index_row,index_col] = float(line.rstrip().split("\t")[1])

        matrix_write(matcov,output["cov"],sorted_samples,sorted_mags)
        # concat the other
        with open(output["table"],"w") as handle_w:
            header = next(open(input["table"][0]))
            handle_w.write(header)
            for file in input["table"]:
                with open(file) as handle:
                    header = next(handle)
                    for line in handle:
                        handle_w.write(line)


rule index_contigs:
    input: "{path}/dmag_spk_scg.fa"
    output: "{path}/dmag_spk_scg.lengths"
    singularity: BWASAMTOOLS_IMG
    shell: """
    samtools faidx {input}
    cut -f1,2 {input}.fai > {output}
    """

rule bedtools:
    input:  bam = "{path}/bam/{sample}_filtered.bam",
            bed = "{path}/dmag_spk_scg.bed",
            genome = "{path}/dmag_spk_scg.lengths",
    output: "{path}/cov/{sample}.orfs.cov"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")
    singularity: BEDTOOLS_IMG
    shell:   "bedtools coverage -a {input.bed} -b {input.bam} -g {input.genome} -mean -sorted> {output}"


rule filter_orf_cov:
    input:   cov = "{path}/cov/{sample}.orfs.cov",
             spk = "{path}/spike_in_def.tsv"
    output:  "{path}/cov/{sample}_filtered.orfs.cov"
    params: derep = SUMMARY,
            prefix = "{path}/cov/{sample}",
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")
    shell : "{SCRIPTS2}/filter_dmags_orf_cov.py {input.cov} {params.derep} {input.spk} {output}"


rule collate_coverage:
    input:   filt = [f"{{path}}/cov/{sample}_filtered.orfs.cov" for sample in R1],
             details = [f"{{path}}/cov/{sample}.orfs.cov" for sample in R1]
    output:  filt = "{path}/drep_95_orfs.cov"
    params: suffix = "_filtered.orfs.cov"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")
    shell : """
    {SCRIPTS}/collate_coverage.py -o {output.filt} -s {params.suffix} -l {input.filt} 
    """

def get_mag_dmags_summary(drep_clu):
    mag_to_dmag = {}
    dmag_to_mags = defaultdict(list)
    with open(drep_clu) as handle:
        _ = next(handle)
        for line in handle:
            sline = line.rstrip().split("\t")
            mag,_,dmag = sline[:3]
            mag_to_dmag[mag] = dmag
            dmag_to_mags[dmag].append(mag)
    return mag_to_dmag,dmag_to_mags

def add_spike_in(mag_to_clu,clu_def,spike_in):
    # currently needed because results summary doesn't contain the spike in 
    spkin = [line.rstrip() for line in open(spike_in)]
    for spk in spkin:
        mag_to_clu[spk] = spk
        clu_def[spk] = [spk]
    return mag_to_clu,clu_def


rule collate_orfs_details:
    input: cov = [f"{{path}}/cov/{sample}.orfs.cov" for sample in R1],
           spk = "{path}/spike_in_def.tsv"
    output: "{path}/drep_95_orfs_detail.cov"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")    
    run:
        mag_to_dmag, dmag_to_mags = get_mag_dmags_summary(SUMMARY)
        mag_to_dmag, dmag_to_mags = add_spike_in(mag_to_dmag, dmag_to_mags, input["spk"])
        with open(output[0],"w") as handle:
            handle.write("%s\n"%"\t".join(["dmag","scg","cov","sample"]))
            for file in input["cov"]:
                dmag_scg_cov = defaultdict(float)
                sample = basename(file).replace(".orfs.cov","")
                for line in open(file):
                    contig,start,end,orf,scg,mag,cov = line.rstrip().split("\t")
                    dmag_scg_cov[(mag_to_dmag[mag],scg)]+=float(cov)
                results = [[dmag,scg,cov] for (dmag,scg),cov in dmag_scg_cov.items()]
                results = sorted(results,key=lambda x:[x[0],x[1]])
                handle.writelines(f"{dmag}\t{scg}\t{cov}\t{sample}\n" for dmag,scg,cov in results)





# --------------- get scg from spike in -----------------------------
rule prodigal:
    input:"{path}.fasta"
    output:
        faa ="{path}.faa",
        fna="{path}.fna",
        gff="{path}.gff"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")
    singularity: PRODIGAL_IMG
    shell: "prodigal-gv -i {input} -a {output.faa} -d {output.fna} -f gff -o {output.gff}"



rule Batch_rpsblast:
    input:   "{path}.faa"
    output:  "{path}.cogs.tsv"
    params:  db = COG_DB
    log:     "{path}_cog.log"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")
    singularity: BLAST_IMG
    shell:   """
             rpsblast -outfmt '6 qseqid sseqid evalue pident length slen qlen' -evalue 0.00001 -query {input} -db {params.db} -out {output} &>log
             """

rule parse_cogs_annotation:
    input:   cog = "{path}.cogs.tsv"
    output:  cog="{path}_cogs_best_hits.tsv"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")                 
    singularity: PYTHONENV_IMG
    shell:   """
             {SCRIPTS}/Filter_Cogs.py {input.cog} --cdd_cog_file {SCG_DATA}/cdd_to_cog.tsv  > {output.cog}
             """

rule extract_SCG_sequences:
    input:  annotation="{filename}_cogs_best_hits.tsv",
            gff="{filename}.gff",
            fna="{filename}.fna"
    output: "{filename}_SCG.fna"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")
    shell:  "{SCRIPTS}/Extract_SCG.py {input.fna} {input.annotation} {SCG_DATA}/scg_cogs_min0.97_max1.03_unique_genera.txt {input.gff}>{output}"


rule bed_orfs:
    input:   gff="{path}.gff"
    output:  bed="{path}.bed"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")
    shell : "{SCRIPTS}/Gff_to_bed.py {input.gff} {output.bed}"




# # ------------------ use sylph ! -------------------
rule create_ref_file:
    output: "{path}/mag_list.txt"
    run:
        with open(output[0],"w") as handle:
            handle.writelines("%s\n"%mag_path for mag_path in DMAGs)


rule sylph_sketch:
    input: "{path}/mag_list.txt"
    output: "{path}/dmag.syldb"
    params: "{path}/dmag"
    threads: 24
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")
    singularity: SYLPH_IMG
    shell: "sylph sketch -c 50 -l {input} -t {threads} -o {params}" 

rule sample_sketch:
    output: "{path}/samples/{sample}.paired.sylsp"
    params: out = "{path}/samples",
            R1 = lambda w:R1[w.sample],
            R2 = lambda w:R2[w.sample],
            sample = "{sample}"
    threads: 4
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")
    singularity: SYLPH_IMG
    shell: "sylph sketch -c 50 -t {threads} -1 {params.R1} -2 {params.R2} -d {params.out} -S {params.sample}"

rule sylph_profile:
    input: db = "{path}/dmag.syldb",
           spl = "{path}/samples/{sample}.paired.sylsp"
    output: "{path}/cov/{sample}.tsv"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")
    threads: 4
    singularity: SYLPH_IMG
    shell: "sylph profile {input.db} {input.spl}  -t {threads} -o {output}"
