#!/bin/bash



####--------------------------------------
##SLURM options
####--------------------------------------
#SBATCH --job-name picard_tools_insert_size_metrics
#SBATCH --account nvastenh_competition_model 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --time 00:30:00
#SBATCH --array=2-63
#SBATCH --chdir=/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter/scripts
#SBATCH --output=/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter/scripts/out/%x_%j_%a.out
#SBATCH --error=/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter/scripts/err/%x_%j_%a.err




####--------------------------------------
## set variables
####--------------------------------------

echo "---1. Setting variables"

workdir=/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter
sample_metadata=$workdir/scripts/array.tsv
sample_name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${sample_metadata} | awk -F'\t' '{print $4}') #SRR
library_type=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${sample_metadata} | awk -F'\t' '{print $10}') #paired or single
target=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${sample_metadata} | awk -F'\t' '{print $6}') #histon mod etc
INPUT_BAM=$workdir/02.mapping/00.Raw/$sample_name.no_MT.sorted.bam
replicate_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${sample_metadata} | awk -F'\t' '{print $9}')
outpath=${workdir}/02.mapping/00.Raw.PicardIS
OUTPUT=${outpath}/$sample_name.picardIS.txt
HIST_OUT=$outpath/$sample_name.picardIS.pdf

# check if directory exist of not, if not create it
echo "Creating directory: $outpath"
mkdir -p "$outpath"

####--------------------------------------
## loading modules
####--------------------------------------

echo "---2. Load modules"
module load picard/3.0.0


####--------------------------------------
## run cutadapt on all fastqs
####--------------------------------------

echo "---3. Running cutadapt on all fastq"


if [[ "$library_type" == "PAIRED" ]]; then

	picard CollectInsertSizeMetrics H=$HIST_OUT I=$INPUT_BAM O=$OUTPUT 

else
	echo "SINGLE END LIBRARY"
	picard CollectInsertSizeMetrics H=$HIST_OUT I=$INPUT_BAM O=$OUTPUT
fi


