# DECODE pipeline

## About
**Long read analysis workflow**

Repository containing [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow for 
- running [Kraken2](https://ccb.jhu.edu/software/kraken/) taxonomy on preprocessed reads
- assembly with [Megahit](https://github.com/voutcn/megahit), [SPAdes](https://github.com/ablab/spades) and [PLASS](https://github.com/soedinglab/plass)
- gene calls with [Prodigal](https://github.com/hyattpd/Prodigal)
- eukaryote classification with [EUKulele](https://github.com/AlexanderLabWHOI/EUKulele)
 
# Setup

## Cloning the repository

```bash
git clone --recurse-submodules <repo https/ssh URL>
```
If you cloned the project but forgot `--recurse-submodules` do
```bash
git submodule update --init --recursive
```
in the cloned repository.


## Sample setup
Provide a `config/sample.tsv` file in the following format

|Sample_ID | sR1 | sR2 | long_reads | group |
|-----------|-----|-----|------------|-------|
| sample1 | DIRECTORY_PATH/sample1_R1.fastq | DIRECTORY_PATH/sample1_R2.fastq | | A |
| sample2 | DIRECTORY_PATH/sample2_R1.fastq | DIRECTORY_PATH/sample2_R2.fastq | | A |
| sample3 | DIRECTORY_PATH/sample3_R1.fastq | DIRECTORY_PATH/sample3_R2.fastq | | B |
| sample4 | DIRECTORY_PATH/sample4_R1.fastq | DIRECTORY_PATH/sample4_R2.fastq | | B |

- `long_reads`: leave column blank if not available
- `group`: sample groupings for co-assembly and co-binning 


## Databases
Downloaded dbs on `2023-09-13` as indicated below 
- from the [Kraken2 database page](https://benlangmead.github.io/aws-indexes/k2)
    - `/hdd0/susbus/databases/kraken2/pluspfp`: [pluspfp](https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20230605.tar.g)
    - `/hdd0/susbus/databases/kraken2/standard_20230913`: [standard](https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230605.tar.gz)

## BRACKEN
- Prior to running, install [bracken](https://github.com/jenniferlu717/Bracken)
- After installation, set *execute* permissions for the `bracken` & `combine_bracken_outputs.py` script
```bash
chmod +x PATH_to_INSTALL_DIR/Bracken/bracken
chmod +x PATH_to_INSTALL_DIR/Bracken/analyses_scripts/combine_bracken_outputs.py
```

## Conda

[Conda user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

```bash
# install miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod u+x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh # follow the instructions
```

Create the main `snakemake` environment

```bash
# create venv
conda env create -f envs/requirements.yaml
conda activate decode``
```


## Configuration

All config files are stored in the folder `config/`

**Important Note(s)**: 
- Edit the paths to 
    - `data_dir`: *eg:* `/ssd0/susbus/socd/data/preprocessed`
    - `results_dir`: *eg:* `/ssd0/susbus/socd/results`
    - `env_dir`: *eg:* `/ssd0/susbus/socd/kraken2/envs`
    - ***`databases`***: *eg:* `/hdd0/susbus/databases/kraken2/pluspfp`

- Provide a `sample.tsv` in the *config* folder. See format above.

All workflows have a `snakemake` profile (`workflow/profiles/slurm/`) including a `snakemake` config file (`config.yaml`) and a `slurm` config file (`slurm.yaml`) for execution on the HPC cluster.

See the `README.md` inside each workflow folder for more information on the executed steps and configuration.


## IMPORTANT
Config files:
- `config/config.yaml`: main config file for all workflows
- `workflow/profiles/slurm/slurm.yaml`: `slurm` config
- `workflow/profiles/slurm/config.yaml`: `snakemake` parameters

Before executing the workflow, check **all** config files listed above, especially lines tagged with `USER_INPUT`.

It is **not** recommended to run the workflow on the access node:
though all computation-intensive steps should be submitted via `slurm`, it is better to avoid doing that especially for (very) long jobs.


## EXECUTION
```bash
# start an interactive session

# activate the main conda env.
conda activate decode

# dry-run
snakemake --profile profiles/slurm --dry-run
# execute w/ slurm
snakemake --profile profiles/slurm
```
