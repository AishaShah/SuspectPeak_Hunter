#!/bin/bash


module load gcc samtools bedtools2 

#---------------------------------------#
# 01. Set up directories and file paths #
#---------------------------------------#
sample_name=$1
bam_file=$2
output_filtered_bam=$3
INTERMEDIATE=$4
threads=$5
mkdir -p ${INTERMEDIATE}
#---------------------------------------#
# 02. Filter                            #
#---------------------------------------#


filter_start_time=$(date +%s)


## get read ids of reads multimapped
echo "---- $(date +'%H:%M:%S')|${sample_name}:01:Getting IDs of Multimapped Reads ------"
multimapped_reads=${INTERMEDIATE}/${sample_name}.ReadIDs_multimapped.out
samtools view "$bam_file" | fgrep "XS:" | cut -f 1 | sort | uniq > ${multimapped_reads}

## remove multimapped reads + reads mapped to repeats
echo "---- $(date +'%H:%M:%S')|${sample_name}:02: Removing Reads overlapping repeats and Multimmaped ------"
samtools view --threads ${threads} -h $bam_file | fgrep -w -vf <(cat ${multimapped_reads}) |  samtools sort --threads ${threads} -n - > ${output_filtered_bam} #${INTERMEDIATE}/$sample_name.no_MT.UM.bam


filter_end_time=$(date +%s)

function calculate_time {
    local start=$1
    local end=$2
    local duration=$((end - start))
    echo "$(date -u -d @${duration} +"%T")"
}

filter_time=$(calculate_time $filter_start_time $filter_end_time)

echo -e "filter"
echo -e "${filter_time}"
