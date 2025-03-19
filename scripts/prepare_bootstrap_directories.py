import pandas as pd
import numpy as np
import argparse
from collections import defaultdict
import os


def load_metadata(file_path):
    return pd.read_csv(file_path, sep='\t')


def create_bootstrap_directories(array_sample_list, round, dir):
    dir_for_round=os.path.join(dir, f"round_{round}","00.mapping")
    os.makedirs(dir_for_round, exist_ok=True)
    for sample in array_sample_list:
    
        # Extract sample names
        original_file = f"02.mapping/00.Raw/{sample}.no_MT.sorted.bam"
        link_name = os.path.join(dir_for_round,f"{sample}.no_MT.sorted.bam")
        os.symlink(os.path.relpath(original_file, start=dir_for_round), link_name)

        #index files
        original_file = f"02.mapping/00.Raw/{sample}.no_MT.sorted.bam.bai"
        link_name = os.path.join(dir_for_round, f"{sample}.no_MT.sorted.bam.bai")
        os.symlink(os.path.relpath(original_file, start=dir_for_round), link_name)
        print(f"Created directory and soft links for round {round}\n Samples: {sample}\n")





def main():
    parser = argparse.ArgumentParser(description='Create bootstrap directories and symlinks for a specific round.')
    parser.add_argument('array_dir', type=str, help='Path to the TSV file containing sample names for the round.')
    parser.add_argument('rounds', type=int, help='Number of Bootstrap Rounds.')
    parser.add_argument('BAM_outdir', type=str, help='Base directory for bootstrap BAM softlinks.')
    args = parser.parse_args()

    for round in range(args.rounds):
        round_id=round+1
        array_file = os.path.join(args.array_dir, f'array_{round_id}.tsv')
        round_samples=pd.read_table(array_file).set_index('sample', drop=False).index.to_list()
        BAM_outdir=args.BAM_outdir
        # Create bootstrap directories and symlinks for the given round
        create_bootstrap_directories(round_samples, round_id, BAM_outdir)

if __name__ == '__main__':
    main()