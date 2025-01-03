#!/bin/bash

# Usage: ./remove_mitochondrial_chromosome.sh <sample> <bam> <intermediate_dir> <output_BAM> <threads>
# Example: ./remove_mitochondrial_chromosome.sh sample_name 8

# Parse command line arguments
SAMPLE=$1
INPUT_BAM=$2
INTERMEDIATE_DIR=$3
BAM_WO_MT=$4
THREADS=$5


# Create intermediate directory
mkdir -p ${INTERMEDIATE_DIR}

# Remove mitochondrial chromosome and sort the BAM file
samtools view --threads ${THREADS} -bS ${INPUT_BAM} -e 'rname != "MT"' > ${INTERMEDIATE_DIR}/${SAMPLE}.no_MT.bam
samtools sort --threads ${THREADS} ${INTERMEDIATE_DIR}/${SAMPLE}.no_MT.bam > ${BAM_WO_MT}

# Index the resulting BAM file
samtools index -@ ${THREADS} ${BAM_WO_MT}

# Remove intermediate files if keep_BAMs_wo_MT is not True
# if [ "${KEEP_BAMS_WO_MT}" != "True" ]; then
#     rm -rf ${INTERMEDIATE_DIR}
# fi


### stats
# Get counts from samtools view -c commands
mapped_reads=$(samtools view -c "${INPUT_BAM}")
removing_MT=$(samtools view -c "${BAM_WO_MT}")
mkdir -p "02.mapping/stats"
# Print the heading and the counts in a single row
echo -e "Sample\tMapped reads\tRemoving MT" > "02.mapping/stats/${SAMPLE}.Removed_MT.bam.stats"
echo -e "${SAMPLE}\t${mapped_reads}\t${removing_MT}" >> "02.mapping/stats/${SAMPLE}.Removed_MT.bam.stats"

echo "Processing of ${SAMPLE} completed."
