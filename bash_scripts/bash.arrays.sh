#!/bin/bash
#SBATCH --account=nvastenh_competition_model
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --job-name=snakemake_array
#SBATCH --array=1,3
#SBATCH --chdir=/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter

# Dynamically set log directory based on the array ID
LOG_DIR="results/config_${SLURM_ARRAY_TASK_ID}/logs"
mkdir -p $LOG_DIR

#SBATCH --output=${LOG_DIR}/snakemake_array_%x_%A_config_%a.out
#SBATCH --error=${LOG_DIR}/snakemake_array_%x_%A_config_%a.err

# Load required modules and activate the conda environment
conda activate /work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/conda_env/trimgalore_env_mamba
module load gcc r-light bedtools2 samtools bowtie2

# Set directories
WORKDIR="/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter"
CONFIG_FILE="$WORKDIR/Config/config_${SLURM_ARRAY_TASK_ID}.yaml"
SNAKEFILE="$WORKDIR/Snakefile"

# Ensure the results directory exists
mkdir -p $RESULTS_DIR



# Run Snakemake for the current array task
snakemake \
  --configfile $CONFIG_FILE \
  --cluster-config cluster.yaml --cluster "sbatch --cpus-per-task={threads} --output={cluster.output} --account={cluster.account} --time={cluster.time}" \
  --jobs 60 \
  --latency-wait 5 \
  --keep-going \
  --rerun-incomplete \
  --snakefile $SNAKEFILE Initialize_SuspectPeak_Hunter


snakemake \
  --configfile $CONFIG_FILE \
  --cluster-config cluster.yaml --cluster "sbatch --cpus-per-task={threads} --output={cluster.output} --account={cluster.account} --time={cluster.time}" \
  --jobs 60 \
  --latency-wait 5 \
  --keep-going \
  --rerun-incomplete \
  --snakefile $SNAKEFILE Validate_SL
