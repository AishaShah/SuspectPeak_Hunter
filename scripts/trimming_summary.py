import os
import sys

def extract_last_line_value(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        last_line = lines[-1].strip()
        if last_line:  # Check if the last line is not empty
            return last_line.split(': ')[-1].strip()
        else:  # If the last line is empty, return the value of the second last line
            second_last_line = lines[-2].strip()
            return second_last_line.split(':\t')[-1].strip()

# Check if directory argument is provided
if len(sys.argv) < 2:
    print("Usage: python script.py <directory>")
    sys.exit(1)

directory = sys.argv[1]

#directory = "/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/01_pipeline/01.Prepare_Input_Data/fastq"

# Define headers
headers = ["Sample", "Total reads processed", "Reads with adapters", "Reads written (passing filters)",
           "Total basepairs processed", "Quality-trimmed", "Total written (filtered)", "Reads Removed"]

print(*headers, sep='\t')

for root, dirs, files in os.walk(directory):
    for file in files:
        if file.endswith("trimming_report.txt"):
            file_path = os.path.join(root, file)
            sample = os.path.splitext(file)[0].replace("_trimming_report", "")
            relevant_stats = []
            with open(file_path, 'r') as f:
                for line in f:
                    if any(metric in line for metric in ['Total reads processed:', 'Reads with adapters:', 'Reads written (passing filters):', 'Total basepairs processed:', 'Quality-trimmed:', 'Total written (filtered):']):
                        relevant_stats.append(line.split(': ')[-1].strip())
            last_line_value = extract_last_line_value(file_path)
            print(sample, *relevant_stats, last_line_value, sep='\t')

