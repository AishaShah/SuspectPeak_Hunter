#!/bin/bash



####--------------------------------------
##SLURM options
####--------------------------------------
#SBATCH --job-name SEACR_thresholds
#SBATCH --account nvastenh_competition_model 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 00:10:00
#SBATCH --array=3-63
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
INPUT_BEDGRAPH=$workdir/03.PeakCalling/00.rawData/DownSampling_3000000/01.genomecov/$sample_name.bedgraph
outpath=${workdir}/03.PeakCalling/00.rawData/DownSampling_3000000/01.genomecov.SEACR.threshold

# check if directory exist of not, if not create it
echo "Creating directory: $outpath"
mkdir -p "$outpath"


# code from: https://github.com/FredHutch/SEACR/blob/master/SEACR_1.3.sh 
awk 'BEGIN{s=1}; {if(s==1){s++}else if(s==2){if($4 > 0){chr=$1; start=$2; stop=$3; max=$4; coord=$1":"$2"-"$3; auc=$4*($3-$2); num=1; s++}}else{if($4 > 0){if(chr==$1 && $2==stop){num++; stop=$3; auc=auc+($4*($3-$2)); if ($4 > max){max=$4; coord=$1":"$2"-"$3}else if($4 == max){split(coord,t,"-"); coord=t[1]"-"$3}}else{print chr"\t"start"\t"stop"\t"auc"\t"max"\t"coord"\t"num; chr=$1; start=$2; stop=$3; max=$4; coord=$1":"$2"-"$3; auc=$4*($3-$2); num=1}}}}' $INPUT_BEDGRAPH > $outpath/$sample_name.auc.bed
cut -f 4,7 $outpath/$sample_name.auc.bed > $outpath/$sample_name.auc

