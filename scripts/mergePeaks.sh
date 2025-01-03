#!/bin/bash
module load gcc bedtools2
# Function to display help message
display_help() {
    echo "Usage: $0 [FILE]..."
    echo "Merge SEACR PeakFiles"
    echo "Example: $0 a.bed b.bed c.bed"
    exit 1
}

# Check if no arguments are provided or if -h or --help is provided
if [[ $# -eq 0 || "$1" == "-h" || "$1" == "--help" ]]; then
    display_help
fi

# Check if mergeBed is available
if ! command -v mergeBed &>/dev/null; then
    echo "Error: mergeBed not found. Please make sure BEDTools is installed and accessible."
    exit 1
fi

# Use awk to append filename, sort, and mergeBed
echo -e "chrom\tstart\tend\tlength\tnum\tnames" 
#awk '{print $1 "\t" $2 "\t" $3 "\t" FILENAME}' "$@" \
awk '{filename=FILENAME; sub(/^.*\//, "", filename); print $1 "\t" $2 "\t" $3 "\t" filename}' "$@" \
| sort -k1,1 -k2,2n \
| mergeBed -i stdin -header -c 4 -o count_distinct,distinct | awk '{print $1 "\t" $2 "\t" $3 "\t" $3-$2 "\t" $4 "\t" $5}' 

