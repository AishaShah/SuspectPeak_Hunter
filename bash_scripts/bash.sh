#!/bin/bash
#SBATCH --account=nvastenh_competition_model
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --job-name=snakemake
#SBATCH --chdir=/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter
#SBATCH --output=logs/snakemake_%x_%j.out
#SBATCH --error=logs/snakemake_%x_%j.err

#conda init bash
conda activate /work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/conda_env/trimgalore_env_mamba
module load gcc r-light bedtools2 samtools bowtie2


# Set the working directory
WORKDIR="/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter"

# Change to the working directory
cd $WORKDIR

snakemake \
  --cluster-config cluster.yaml Initialize_SuspectPeak_Hunter --cluster "sbatch --cpus-per-task={threads} --output={cluster.output} --account={cluster.account} --time={cluster.time}" \
  --jobs 60 \
  --latency-wait 5 \
  --keep-going \
  --rerun-incomplete \
  --snakefile Snakefile

snakemake \
  --cluster-config cluster.yaml Validate_SL --cluster "sbatch --cpus-per-task={threads} --output={cluster.output} --account={cluster.account} --time={cluster.time}" \
  --jobs 60 \
  --latency-wait 5 \
  --keep-going \
  --rerun-incomplete \
  --snakefile Snakefile
