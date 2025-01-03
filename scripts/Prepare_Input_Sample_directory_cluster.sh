#!/bin/bash

# Check if the correct number of arguments is provided
if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <path_to_array_tsv> <target_directory_fastq> <target_directory_trimmed> <target_directory_bam>"
    exit 1
fi

# Get the arguments
array_tsv=$1
target_dir_fq=$2
target_dir_trimmed_fq=$3
target_dir_bam=$4

# Create the directory if it doesn't exist
mkdir -p $target_dir_fq
mkdir -p $target_dir_trimmed_fq
mkdir -p $target_dir_bam

# Initialize a flag to track if all files are successfully copied
all_files_copied=true

# Function to check if a file is a symbolic link and is accessible
check_path() {
    local file_path=$1
    local file_name=$2
    local file_type=$3
    if [[ ! -e "$file_path" ]]; then
        echo "$file_type:Error: $file_path is not accessible. Recopy Please!"
        all_files_copied=false
    else
        echo "$file_type:$file_path is copied successfully!"
    fi
}


copy_file() {
    local path_fq=$1
    local file_path=$2
    local file_name=$3
    local file_type=$4
    if [[ -e "$file_path" ]]; then
        cp "${path_fq}/${file_name}" "${file_path}"
    else
        echo "$file_type:Skipping:$file_path already exists"
    fi
}



# Read the sample sheet and create soft links
while IFS=$'\t' read -r projectID ctrl shtname sample Stage Target Target_Group control rep type path_fq path_fq_trim R1 R2 path_bam Fwd_adaptor Rev_adaptor
do
    # Skip the header line
    if [[ "$projectID" == "projectID" ]]; then
        continue
    fi
    # Create soft links based on file type

    # fastq file
    if [[ -n "$path_fq" && "$path_fq" != "NA" ]]; then
        
        # Echo the commands for paired-end reads

        if [[ "$type" == "PAIRED" ]]; then
            #copy_file "${path_fq}/${R1}" "${target_dir_fq}/${sample}_1.fastq.gz" "${R1}" "PE"
            #copy_file "${path_fq}/${R2}" "${target_dir_fq}/${sample}_2.fastq.gz" "${R2}" "PE"
            cp "${path_fq}/${R1}" "${target_dir_fq}/${sample}_1.fastq.gz"
            cp "${path_fq}/${R2}" "${target_dir_fq}/${sample}_2.fastq.gz"
            check_path "${target_dir_fq}/${sample}_1.fastq.gz" "${R1}" "PE"
            check_path "${target_dir_fq}/${sample}_2.fastq.gz" "${R2}" "PE"
        # Echo the command for single-end read
        elif [[ "$type" == "SINGLE" ]]; then
            cp "${path_fq}/${R1}" "${target_dir_fq}/${sample}.fastq.gz"
            #copy_file "${path_fq}/${R1}" "${target_dir_fq}/${sample}.fastq.gz" "${R1}" "SE"
            check_path "${target_dir_fq}/${sample}.fastq.gz" "${R1}" "SE"
        fi
    fi

    # fastq trimmed files
    if [[ -n "$path_fq_trim" && "$path_fq_trim" != "NA" ]]; then
        
        # Echo the commands for paired-end reads
        if [[ "$type" == "PAIRED" ]]; then

            cp "${path_fq_trim}/${R1}" "${target_dir_trimmed_fq}/${R1}"
            cp "${path_fq_trim}/${R2}" "${target_dir_trimmed_fq}/${R2}"
            check_path "${target_dir_trimmed_fq}/${R1}" "${R1}" "PE"
            check_path "${target_dir_trimmed_fq}/${R2}" "${R2}" "PE"
        # Echo the command for single-end read
        elif [[ "$type" == "SINGLE" ]]; then

            cp "${path_fq_trim}/${R1}" "${target_dir_trimmed_fq}/${R1}"
            check_path "${target_dir_trimmed_fq}/${R1}" "${R1}" "SE"
        fi
    fi    

    # BAM files
    if [[ -n "$path_bam" && "$path_bam" != "NA" ]]; then
        echo "Creating soft links for BAM files"
        cp "$path_bam/${sample}.bam" "$target_dir_bam/${sample}.bam"
        check_path "${target_dir_bam}/${sample}.bam" "${sample}" "BAM"
    fi

done < <(cat $array_tsv)

# Report final status
if $all_files_copied; then
    echo "All files have been successfully copied and are accessible." 
else
    echo "There were errors linking some files. Please check the log above."
fi




