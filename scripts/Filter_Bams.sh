#!/bin/bash


module load gcc samtools bedtools2 

#---------------------------------------#
# 01. Set up directories and file paths #
#---------------------------------------#
sample_name=$1
sample_type=$2
bam_file=$3
regions_comp_to_repeat=$4
genome=$5
output_filtered_bam=$6
INTERMEDIATE=$7
threads=$8
mkdir -p ${INTERMEDIATE}
#---------------------------------------#
# 02. Filter                            #
#---------------------------------------#


filter_start_time=$(date +%s)


echo "---- $(date +'%H:%M:%S')|${sample_name}:00:Getting IDs of Reads having less than 30bp out of repeat region ------"
repeat_reads=${INTERMEDIATE}/${sample_name}.ReadIDs_overlapping_Repeats.out

## get read ids of reads overlapping repeats
bedtools intersect -a $bam_file -b "$regions_comp_to_repeat" -bed -wao|\
  awk '{if($16<30) print $4}' |\
  sed 's/\/[12]$//' |\
  sort |\
  uniq > ${repeat_reads}


## get read ids of reads multimapped
echo "---- $(date +'%H:%M:%S')|${sample_name}:01:Getting IDs of Multimapped Reads ------"
multimapped_reads=${INTERMEDIATE}/${sample_name}.ReadIDs_multimapped.out
samtools view "$bam_file" | fgrep "XS:" | cut -f 1 | sort | uniq > ${multimapped_reads}

## remove multimapped reads + reads mapped to repeats
echo "---- $(date +'%H:%M:%S')|${sample_name}:02: Removing Reads overlapping repeats and Multimmaped ------"
samtools view -h $bam_file | fgrep -w -vf <(cat ${multimapped_reads} ${repeat_reads}) | samtools view -Sb - > ${INTERMEDIATE}/$sample_name.UM.no_repeats.bam

## fixmate and sorting 

if [[ "$sample_type" == "PAIRED" ]]; then
  # echo "---- $(date +'%H:%M:%S')|${sample_name}:03:fixmate and sorting ------"
  samtools view -h ${INTERMEDIATE}/$sample_name.UM.no_repeats.bam |  samtools sort --threads ${threads} -n - | samtools fixmate - ${INTERMEDIATE}/$sample_name.UM.no_repeats.fixmate.bam
  samtools view -h -f 1 ${INTERMEDIATE}/$sample_name.UM.no_repeats.fixmate.bam > ${output_filtered_bam}
  
  fixmate=$(samtools view -c "${INTERMEDIATE}/${sample_name}.UM.no_repeats.fixmate.bam")

elif [[ "$sample_type" == "SINGLE" ]]; then
  samtools view -h ${INTERMEDIATE}/$sample_name.UM.no_repeats.bam |  samtools sort --threads ${threads} - > ${output_filtered_bam} 
  fixmate="NA"

fi


### stats
# Get counts from samtools view -c commands
mapped_reads=$(samtools view -c "$bam_file")
removing_MM_and_ROR=$(samtools view -c "${INTERMEDIATE}/${sample_name}.UM.no_repeats.bam")
properly_paired=$(samtools view -c "${output_filtered_bam}")

mkdir -p "02.mapping/stats"
# Print the heading and the counts in a single row
echo -e "Sample\tMapped reads\tRemoving MM and ROR\tFixmate\tProperly Paired" > "02.mapping/stats/${sample_name}.Filtered.bam.stats"
echo -e "${sample_name}\t${mapped_reads}\t${removing_MM_and_ROR}\t${fixmate}\t${properly_paired}" >> "02.mapping/stats/${sample_name}.Filtered.bam.stats"


## remove multimapped reads + reads mapped to repeats
# echo "---- $(date +'%H:%M:%S')|${sample_name}:02: Removing Reads overlapping repeats and Multimmaped ------"
#samtools view -h $bam_file | fgrep -w -vf <(cat ${multimapped_reads} ${repeat_reads}) | samtools sort --threads ${threads} -n - | samtools fixmate - ${INTERMEDIATE}/$sample_name.UM.no_repeats.bam
## remove reads that are over-softclipped 
# echo "---- $(date +'%H:%M:%S')|${sample_name}:03:Removing reads having >30bp softclipped at any end ------"
# samtools view -h ${INTERMEDIATE}/$sample_name.UM.no_repeats.bam | samclip --ref ${genome} --max 30 --progress 0 | samtools sort -n - | samtools fixmate - ${INTERMEDIATE}/$sample_name.UM.no_repeats.SC30.fixmate.bam
# samtools view -h -f 1 ${INTERMEDIATE}/$sample_name.UM.no_repeats.SC30.fixmate.bam > ${output_filtered_bam}

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
