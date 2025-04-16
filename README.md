This Snakemake workflow takes paired-end whole-genome bisulfite sequencing (WGBS) data and processes it using
BISulfite-seq CUI Toolkit (BISCUIT).

BISCUIT was written to perform alignment, DNA methylation and mutation calling, and allele specific methylation from
bisulfite sequencing data (https://huishenlab.github.io/biscuit/).

Download BISCUIT here: https://github.com/huishenlab/biscuit/releases/latest.

# Components of the workflow
The following components are generally in order, but may run in a different order, depending on exact dependencies
needed.
  + [default off] Generate asset files used during QC related rules
  + [default off] Modify and index reference genome to include methylation controls (lambda phage and pUC19)
  + [default off] Trim FASTQ files
  + [default off] Run Fastq Screen in bisulfite mode
  + Run FastQC on raw FASTQ files
  + Alignment, duplicate marking, and indexing of input data (biscuitSifter pipeline)
  + Samtools flagstat of input data
  + Methylation information extraction (BED Format)
  + Merge C and G beta values in CpG dinucleotide context
  + [default off] SNP and epiBED extraction
  + [default off] Run Preseq on aligned BAM
  + MultiQC with BICUIT QC modules specifically for methylation data
  + [default off] Generate plots of the observed / expected coverage ratio for different genomic features
  + [default off] Generate percentage of covered CpGs and CpG island coverage figures
  + [default off] Find coverage uniformity across genome
  + [default off] Plot percentage of genome covered
  + [default off] Find average methylation values in bins across genome
  + [default off] Find average methylation values in bins centered on specified regions
  + [default off] QC methylated and unmethylated controls
  + [default off] Generate bigWigs from BISCUIT BED files

Many options can be easily specified in the `config.yaml`! Otherwise, the commands in the Snakefile can also be modified
to meet different needs.

# Dependencies

The following dependencies are downloaded when running with `--use-conda`, otherwise you must have these in your PATH.
| Package        | Conda Version Downloaded | Notes |
|:---------------|:------------------------:|:------|
| `snakemake`    | 9.3+                     | Needed before running pipeline |
| `biscuit`      | 1.6.1                    |       |
| `htslib`       | 1.21                     |       |
| `samtools`     | 1.21                     |       |
| `dupsifter`    | 1.3.0                    |       |
| `parallel`     | 20230322                 |       |
| `bedtools`     | 2.30.1                   |       |
| `perl`         | 5.32.1                   | Only required if running `beta_bigwigs` |
| `ucsc-bedgraphtobigwig` | 455             | Only required if running `beta_bigwigs` |
| `preseq`       | 3.2.0                    | Must be compiled with htslib enabled |
| `fastqc`       | 0.12.1                   |       |
| `trim_galore`  | 0.6.10                   |       |
| `fastq_screen` | 0.16.0                   | Only required if running `fastq_screen` |
| `bismark`      | 0.24.2                   | Only required if running `fastq_screen` |
| `pigz`         | 2.8                      |       |
| `python`       | 3.13.3                   |       |
| `pandas`       | 2.2.3                    |       |
| `numpy`        | 2.2.4                    |       |
| `matplotlib`   | 3.10.1                   |       |
| `seaborn`      | 0.13.2                   |       |
| `multiqc`      | 1.28                     |       |
| `R`            | 4.4.3                    |       |
| `tidyverse`    | 2.0.0                    | Only required for plotting methylation controls |
| `ggplot2`      | 3.5.2                    | Only required for plotting methylation controls |
| `patchwork`    | 1.3.0                    | Only required for plotting methylation controls |
| `viridislite`  | 0.4.2                    | Only required for plotting methylation controls |

Two things of note, 1) it is easiest when working with `snakemake` to install `mamba` using `conda` when running with
`--use-conda`, and 2) it is preferable to install `snakemake` using `conda`, rather than using a module. This is due to
potential conflicts between packages (such as `matplotlib`) that can be found in the snakemake module's python
distrubtion and the conda installed python distribution.

# Running the workflow
For ease of reference, the configuration file `config/config.yaml` will be referred to throughout as the file to define
any configuration needed for your pipeline run. That said, you can copy this config file to another file and use that
config file in your pipeline with `snakemake --configfile /my/new/config.yaml` or by changing the `CONFIG_FILE` variable
in the SLURM submit script.

+ Clone the repo (https://github.com/huishenlab/Biscuit_Snakemake_Workflow/tree/master).
  + SSH: `git clone git@github.com:huishenlab/Biscuit_Snakemake_Workflow.git`
  + HTTPS: `git clone https://github.com/huishenlab/Biscuit_Snakemake_Workflow.git`

+ Place *gzipped* FASTQ files into `raw_data/`. Alternatively, you can specify the location of your *gzipped* FASTQ
files in `config/config.yaml`.

+ Replace the example `config/samples.tsv` with your own sample sheet containing:
  + A row for each sample
  + The following three columns for each row:
    + A. `sample`
    + B. `fq1` (name of R1 file for `sample` in your raw data directory)
    + C. `fq2` (name of R2 file for `sample` in your raw data directory)
    + D. Any other columns included are ignored
  + Note, you can either edit `config/samples.tsv` in place or specify the path to your sample sheet in
  `config/config.yaml`. If you create your own sample sheet, make sure to include the header line as is seen in the
  example file.

+ Modify `config/config.yaml` to specify the appropriate
  + Reference genome
  + Biscuit index
  + Biscuit QC assets (https://github.com/huishenlab/biscuit/releases/latest)
  + Toggle optional workflow components
  + Set other run parameters in `config/config.yaml`
  + Turn on optional rules in `config/config.yaml` (change from False to True)
  + If you are using environmental modules on your system, you can set the locations in the corresponding location. By
  default, the pipeline will use `conda`/`mamba` to download the required packages. Note, if using the modules and a
  module is not available, snakemake gives a warning but will run successfully *as long as the required executables are
  in the path*.

+ Modify SLURM submit script as needed (new config file in `CONFIG_FILE`, etc.).

+ Then submit the workflow to an HPC using something similar to `bin/run_snakemake_workflow.slurm` (e.g.,
`sbatch bin/run_snakemake_workflow.slurm`). `bin/run_snakemake_workflow.slurm` works for a SLURM queue
system. A PBS/Torque version is available in a previous release on GitHub for those who need it.

# After the workflow

+ The output files in `config['output_directory']/analysis/pileup/` may be imported into a `BSseq` object using
`bicuiteer::readBiscuit()`.
+ `config['output_directory']/analysis/multiqc/multiqc_report.html` contains the methylation-specific BISCUIT QC modules
(https://huishenlab.github.io/biscuit/docs/alignment/QC.html)

# Test dataset

This workflow comes with a working example dataset. To test the smakemake workflow on your system, place the 10
FASTQ files in `bin/working_example_dataset` into `raw_data/` and use the default `config/samples.tsv` sample sheet.
These example files can be mapped to the human genome.

# Example default workflow - 1 sample
![workflow diagram](bin/DAGs/dag.pdf)

# Helpful snakemake commands for debugging a workflow
For more information on Snakemake: https://snakemake.readthedocs.io/en/stable/

  + Do a test run: `snakemake -npr`
  + Unlock directory after a manually aborted run: `snakemake --unlock --cores 1`
  + Create a workflow diagram for your run: `snakemake --dag | dot -Tpng > my_dag.png`
  + Run pipeline from the command line: `snakemake --use-conda --cores 1`
