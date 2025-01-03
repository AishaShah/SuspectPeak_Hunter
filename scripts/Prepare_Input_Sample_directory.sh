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

# Initialize a flag to track if all files are successfully linked
all_files_linked=true

# Function to check if a file is a symbolic link and is accessible
check_symlink() {
    local file_path=$1
    local file_name=$2
    local file_type=$3
    if [[ ! -L "$file_path" || ! -e "$file_path" ]]; then
        echo "$file_type:Error:Soft link for $file_path is not accessible."
        all_files_linked=false
    else
        echo "$file_type:Soft link for $file_path is correct and accessible."
    fi
}

# Read the sample sheet and create soft links
while IFS=$'\t' read -r projectID ctrl shtname sample control rep type path_fq path_fq_trim R1 R2 path_bam
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
            ln -sr "${path_fq}/${R1}" "${target_dir_fq}/${R1}"
            ln -sr "${path_fq}/${R2}" "${target_dir_fq}/${R2}"
            check_symlink "${target_dir_fq}/${R1}" "${R1}" "PE"
            check_symlink "${target_dir_fq}/${R2}" "${R2}" "PE"
        # Echo the command for single-end read
        elif [[ "$type" == "SINGLE" ]]; then
            ln -sr "${path_fq}/${R1}" "${target_dir_fq}/${R1}" 
            check_symlink "${target_dir_fq}/${R1}" "${R1}" "PE"
        fi
    fi

    # fastq trimmed files
    if [[ -n "$path_fq_trim" && "$path_fq_trim" != "NA" ]]; then
        
        # Echo the commands for paired-end reads
        if [[ "$type" == "PAIRED" ]]; then

            ln -sr "${path_fq_trim}/${R1}" "${target_dir_trimmed_fq}/${R1}"
            ln -sr "${path_fq_trim}/${R2}" "${target_dir_trimmed_fq}/${R2}"
            check_symlink "${target_dir_trimmed_fq}/${R1}" "${R1}" "PE"
            check_symlink "${target_dir_trimmed_fq}/${R2}" "${R2}" "PE"
        # Echo the command for single-end read
        elif [[ "$type" == "SINGLE" ]]; then

            ln -sr "${path_fq_trim}/${R1}" "${target_dir_trimmed_fq}/${R1}"
            check_symlink "${target_dir_trimmed_fq}/${R1}" "${R1}" "SE"
        fi
    fi    

    # BAM files
    if [[ -n "$path_bam" && "$path_bam" != "NA" ]]; then
        echo "Creating soft links for BAM files"
        ln -sr "$path_bam/${sample}.bam" "$target_dir_bam/${sample}.bam"
        check_symlink "${target_dir_bam}/${sample}.bam" "${sample}" "BAM"
    fi

done < <(cat $array_tsv)

# Report final status
if $all_files_linked; then
    echo "All files have been successfully linked and are accessible." 
else
    echo "There were errors linking some files. Please check the log above." 
fi
