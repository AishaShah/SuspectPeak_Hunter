#!/bin/bash

# rename this script to bam2bed
module load gcc samtools bedtools2 

#---------------------------------------#
# 01. Set up directories and file paths #
#---------------------------------------#

sample_name=$1
bam_file=$2
intermediate=$3
bamtobed_output=$4
threads=$5
sample_type=$6

mkdir -p $intermediate
#----------------------------------------------#
# 02.generate bed file using bedtools bamtobed #
#----------------------------------------------#


if [[ "$sample_type" == "PAIRED" ]]; then
    echo "---- ${sample_name}:${sample_type}:00:sorting bam and converting to bed ------"
    samtools sort -n --threads ${threads} -T ${intermediate}/${sample_name} ${bam_file} | bedtools bamtobed -bedpe > ${intermediate}/${sample_name}.bamtobed
    echo "---- ${sample_name}:${sample_type}:01:removeing reads having a large insert size i.e read pairs spanning greater than 1000 bp  ------"
    # remove reads having a large insert size i.e read pairs spanning greater than 1000 bp
    awk '$1==$4 && $6-$2 < 1000 {print $0}' ${intermediate}/${sample_name}.bamtobed > ${intermediate}/${sample_name}.clean.bamtobed
    echo "---- ${sample_name}:${sample_type}:02:merging regions spanned by both reads of a pair  ------"
    # merge regions spanned by both reads of a pair 
    cut -f 1,2,6 ${intermediate}/${sample_name}.clean.bamtobed | sort -k1,1 -k2,2n -k3,3n > ${bamtobed_output}

else
    echo "---- ${sample_name}:${sample_type}:00:sorting bam and converting to bed ------"
    bedtools bamtobed -i ${bam_file} | sort -k1,1 -k2,2n -k3,3n | awk '{print $1 "\t" $2 "\t" $3}' > ${bamtobed_output}
fi