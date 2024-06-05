#!/bin/bash


module load gcc samtools bedtools2 

#---------------------------------------#
# 01. Set up directories and file paths #
#---------------------------------------#

sample_name=$1
bam_file=$2
intermediate=$3
bamtobed_output=$4
threads=$5

mkdir -p $intermediate
#----------------------------------------------#
# 02.generate bed file using bedtools bamtobed #
#----------------------------------------------#

echo "---- ${sample_name}:00:sorting bam and converting to bed ------"
samtools sort -n --threads ${threads} -T ${intermediate}/${sample_name} ${bam_file} | bedtools bamtobed -bedpe > ${intermediate}/${sample_name}.bamtobed
echo "---- ${sample_name}:01:removeing reads having a large insert size i.e read pairs spanning greater than 1000 bp  ------"
# remove reads having a large insert size i.e read pairs spanning greater than 1000 bp
awk '$1==$4 && $6-$2 < 1000 {print $0}' ${intermediate}/${sample_name}.bamtobed > ${intermediate}/${sample_name}.clean.bamtobed
echo "---- ${sample_name}:02:merging regions spanned by both reads of a pair  ------"
# merge regions spanned by both reads of a pair 
cut -f 1,2,6 ${intermediate}/${sample_name}.clean.bamtobed | sort -k1,1 -k2,2n -k3,3n > ${bamtobed_output}

