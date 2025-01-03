#!/bin/bash

####--------------------------------------
## SLURM options
####--------------------------------------
#SBATCH --job-name=pheatmaps
#SBATCH --account=nvastenh_competition_model 
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem 100G
#SBATCH --chdir=/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter/scripts
#SBATCH --output=/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter/scripts/out/%x_%j_%a.out
#SBATCH --error=/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter/scripts/err/%x_%j_%a.err

# Load the R module if necessary (module names might differ based on your system)
module load r

# Run the R script
export R_LIBS_USER=/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/R_libs
mkdir -p "/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/R_libs"
Rscript /work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter/scripts/pheatmaps.R

