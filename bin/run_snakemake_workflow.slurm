#!/bin/bash
#SBATCH --export=NONE
#SBATCH -J SNAKEMASTER
#SBATCH -o logs/workflows/workflow_output-%j.log
#SBATCH --ntasks 1
#SBATCH --time 120:00:00
#SBATCH --mem=8G

CONFIG_FILE=config/config.yaml

cd ${SLURM_SUBMIT_DIR}

mkdir -p logs/workflows

# Add snakemake to PATH here
if [[ `which snakemake 2>&1 /dev/null` ]]; then
    snakemake_module="bbc2/snakemake/snakemake-8.25.2"

    module load $snakemake_module
fi

# save DAG job file with time stamp
TIME=$(date "+%Y-%m-%d_%H.%M.%S")
snakemake --configfile ${CONFIG_FILE} --dry-run         > logs/workflows/workflow_${TIME}.txt
snakemake --configfile ${CONFIG_FILE} --dag | dot -Tpng > logs/workflows/workflow_${TIME}.png

# Default to using conda, if using environment modules, then replace --use-conda with --use-envmodules
# Note, this requires downloading mamba (conda install -n base -c conda-forge mamba)
snakemake \
    --printshellcmds \
    --latency-wait 20 \
    --keep-going \
    --use-conda \
    --jobs 20 \
    --configfile ${CONFIG_FILE} \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "mkdir -p logs/{rule}; sbatch \
        --export=ALL \
        --nodes 1 \
        --ntasks-per-node {threads} \
        --mem={resources.mem_gb}G \
        -t {resources.time} \
        -o logs/{rule}/{rule}-%j.log"
