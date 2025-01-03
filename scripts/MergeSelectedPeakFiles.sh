#!/bin/bash
module load gcc bedtools2

# Function to display help message
display_help() {
    echo "Usage: $0 [SAMPLE_FILE] [BED_FILE_DIR] [BED_FILES_SUFFFIX] [merged_regions_binary_data] [merged_regions_collapsed]"
    echo "Merge SEACR PeakFiles for samples listed in SAMPLE_FILE"
    echo "Example: $0 sample_info.txt 00.peaks/00.raw"
    exit 1
}

# Check if correct number of arguments are provided
if [[ $# -ne 5 || "$1" == "-h" || "$1" == "--help" ]]; then
    display_help
fi

SAMPLE_FILE=$1
BED_FILE_DIR=$2
BED_FILES_SUFFIX=$3
merged_regions_binary_data=$4
merged_regions_collapsed=$5
# Check if sample file exists
if [[ ! -f "$SAMPLE_FILE" ]]; then
    echo "Error: Sample file $SAMPLE_FILE not found."
    exit 1
fi

# Check if bed file directory exists
if [[ ! -d "$BED_FILE_DIR" ]]; then
    echo "Error: BED file directory $BED_FILE_DIR not found."
    exit 1
fi

# Check if mergeBed is available
if ! command -v mergeBed &>/dev/null; then
    echo "Error: mergeBed not found. Please make sure BEDTools is installed and accessible."
    exit 1
fi

# Extract sample names from the sample file
SAMPLES=$(awk 'NR > 1 {print $4 "'.$BED_FILES_SUFFIX'"}' "$SAMPLE_FILE" | sort | uniq)

# Initialize a variable to store matched BED paths

BED_FILES=$(for SAMPLE in $SAMPLES; do echo "$BED_FILE_DIR/$SAMPLE"; done)
# Check if each sample's BED file exists in the specified directory

for BED_FILE in $BED_FILES; do
    if [[ ! -f "$BED_FILE" ]]; then
        echo "Warning: BED file not found: $BED_FILES"
    fi
done



# binary
bedtools multiinter -i $BED_FILES -header | \
{
    # Read the first line as the header
    read header

    # Print the header
    echo "$header"

    # Sort the remaining lines
    sort -k1,1V -k2,2n
}  > $merged_regions_binary_data


# Use awk to append filename, sort, and mergeBed
awk '{filename=FILENAME; sub(/^.*\//, "", filename); print $1 "\t" $2 "\t" $3 "\t" filename}' $BED_FILES | \
sort -k1,1V -k2,2n | \
mergeBed -i stdin -header -c 4 -o count_distinct,distinct | \
awk '{print $1 "\t" $2 "\t" $3 "\t" $3-$2 "\t" $4 "\t" $5}' > $merged_regions_collapsed

: '
# JUST TO GET Quick SUMMARY
for i in {1..10}; do
    echo "Results for array_$i:"
    scripts/mergeSelectedPeakFiles.sh "04.bootstrapping_arrays/array_$i.tsv" "03.PeakCalling/00.rawData/DownSampling_3000000/02.peaks" "stringent.bed" | \
    awk {print $5} | sort -k1,1n | uniq -c
    echo ""  # Add an empty line for better readability
done


for i in {1..10}; do
    echo "Results for array_$i:"
    awk {print $5} round_$i.rawbam.mergedPeaks.tab | sort -k1,1n | uniq -c
    echo ""  # Add an empty line for better readability
done
'