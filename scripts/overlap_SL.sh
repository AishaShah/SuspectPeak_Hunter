#!/bin/bash

# Ensure the script exits if any command fails and enable debugging
set -e

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <directory_path> <output_file>"
    exit 1
fi

# Get the input parameters
dir_path="$1"
output_file="$2"

# Ensure the provided directory path exists
if [ ! -d "$dir_path" ]; then
    echo "Error: Directory $dir_path does not exist."
    exit 1
fi

# Process the files and merge intervals, then save to the output file
awk '
{
    # Extract the filename from the full path
    filename = FILENAME
    # Remove the directory path
    sub(/^.*\//, "", filename)
    # Remove the ".SL1.bed" suffix
    sub(/\.SL1\.bed$/, "", filename)
    # Print the desired fields along with the modified filename
    print $1 "\t" $2 "\t" $3 "\t" $4 "\t" filename
}
' "$dir_path"/*.bed |
# Sort the output alphanumerically by the first column and numerically by the second column
sort -k1,1V -k2,2n |
# Use mergeBed to merge overlapping intervals, count distinct filenames, and add a header
mergeBed -i stdin -header -c 5,5 -o count_distinct,distinct > "$output_file"

echo "Processed data saved to $output_file"
