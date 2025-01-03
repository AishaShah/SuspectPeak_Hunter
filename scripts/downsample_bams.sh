#!/bin/bash

# Load necessary modules
module load gcc samtools

# Function to display the help message
display_help() {
    echo "Usage: $0 <BAM_FILE> <Number_of_Reads for PE> <Number_of_Reads for SE> <library_type> <Output_File> <threads>"
    echo
    echo "Description:"
    echo "  Downsamples a BAM file to a specified number of reads."
    echo "  This script calculates the fraction of reads to retain in order to achieve the desired number of reads,"
    echo "  and then generates a downsampled BAM file using that fraction."
    echo
    echo "Example:"
    echo "  $0 Input/a.bam 10000 Downsample/a.bam 8"
    exit 1
}

# Check if no arguments are provided or if -h or --help is provided
if [[ $# -eq 0 || "$1" == "-h" || "$1" == "--help" ]]; then
    display_help
fi

# Assign command-line arguments to variables
bam_file=$1
desired_reads_PE=$2
desired_reads_SE=$3
sample_type=$4
output_file=$5
threads=$6
# Calculate the fraction of reads to retain
if [[ "$sample_type" == "PAIRED" ]]; then
    desired_reads=$desired_reads_PE
    echo "$bam_file:$sample_type:$desired_reads"
else
    desired_reads=$desired_reads_SE
    echo "$bam_file:$sample_type:$desired_reads"
fi



# Calculate the total number of reads in the BAM file
total_reads=$(samtools idxstats "$bam_file" | cut -f3 | awk '{ total += $1 } END { print total }')

# Check if downsampling is needed
if (( total_reads > desired_reads )); then
    # Calculate the fraction for downsampling
    fraction=$(awk -v target_reads="$desired_reads" -v total_reads="$total_reads" 'BEGIN { print target_reads / total_reads }')
    echo "Downsampling with fraction: $fraction"
    
    # Downsample the BAM file
    samtools view --threads "$threads" -b -s "$fraction" "$bam_file" > "$output_file"
    echo "Downsampling complete. Output written to $output_file"
else
    # Copy the original BAM file to the output
    cp "$bam_file" "$output_file"
    echo "Number of reads is less than or equal to desired_reads. Copying original BAM to $output_file"
fi