import pandas as pd
import numpy as np
import argparse
from collections import defaultdict
import os


def load_metadata(file_path):
    return pd.read_csv(file_path, sep='\t')

def is_valid_selection(sample, selected_targets):
    target_stage = (sample['Target'], sample['Stage'], sample['projectID'])
    if any(sample['Target'] == first for (first, _, _) in selected_targets):
        return False
    return target_stage not in selected_targets

def select_samples_for_target_group(group_samples, num_samples, selected_targets):
    selected_samples = []
    attempts = 0
    while len(selected_samples) < num_samples and attempts < 100:
        sample = group_samples.sample(1).iloc[0]
        if is_valid_selection(sample, selected_targets):
            selected_samples.append(sample)
            target_stage = (sample['Target'], sample['Stage'], sample['projectID'])
            selected_targets.add(target_stage)
        attempts += 1
    return selected_samples

def bootstrap_samples(metadata, rounds, samples_per_group):
    results = []
    target_groups = metadata['Target_Group'].unique()

    for _ in range(rounds):
        selected_targets = set()
        round_samples = []
        
        for target_group in target_groups:
            group_samples = metadata[metadata['Target_Group'] == target_group]
            selected_samples = select_samples_for_target_group(group_samples, samples_per_group, selected_targets)
            round_samples.extend(selected_samples)
        
        results.append(pd.DataFrame(round_samples))
    
    return results

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

def check_and_generate_arrays(args, metadata):
    """
    Check for existing array files and generate missing ones if necessary.

    Parameters:
        args (Namespace): Parsed command-line arguments.
        metadata (DataFrame): Loaded metadata for bootstrapping.
    """
    # Expected array file paths
    expected_files = [
        os.path.join(args.array_dir, f'array_{round_id}.tsv')
        for round_id in range(1, args.rounds + 1)
    ]
    
    # Check which files exist
    existing_files = [file for file in expected_files if os.path.exists(file)]
    missing_files = [file for file in expected_files if not os.path.exists(file)]

    if len(existing_files) == len(expected_files):
        print("All array files already exist:")
        for file in existing_files:
            print(f"  - {file}")
        print("No regeneration needed.")
        return

    if missing_files:
        print("The following array files are missing and will be regenerated:")
        for file in missing_files:
            print(f"  - {file}")
    
    print("Regenerating arrays...")
    results = bootstrap_samples(metadata, args.rounds, args.samples_per_group)
    
    # Save regenerated arrays
    for i, round_samples in enumerate(results):
        round_id = i + 1
        output_file = os.path.join(args.array_dir, f'array_{round_id}.tsv')
        round_samples.to_csv(output_file, sep='\t', index=False)
        print(f'Saved Round {round_id} samples to {output_file}')
    
        #sample_names = list(pd.read_table(output_file).set_index('sample', drop=False).index)
        #create_bootstrap_directories(dir=args.bootstrap_dir,array_sample_list=sample_names, round=round_id)


def main():
    parser = argparse.ArgumentParser(description='Perform bootstrapping on sample metadata.')
    parser.add_argument('metadata_file', type=str, help='Path to the metadata file.')
    parser.add_argument('rounds', type=int, help='Number of rounds for bootstrapping.')
    parser.add_argument('samples_per_group', type=int, help='Number of samples to select from each Target_Group in each round.')
    parser.add_argument('array_dir', type=str, help='Directory to save the output TSV files.')
    parser.add_argument('--seed', type=int, default=None, help='Seed for random number generation.')
    #parser.add_argument('bootstrap_dir', type=str, default=None, help='Directory for BAM file shortcuts of each bootstrap.')

    args = parser.parse_args()
    # Set the random seed for reproducibility
    if args.seed is not None:
        np.random.seed(args.seed)
   
    # Ensure the output directory exists
    os.makedirs(args.array_dir, exist_ok=True)
    
    metadata = load_metadata(args.metadata_file)
    check_and_generate_arrays(args, metadata)

        
        

if __name__ == '__main__':
    main()
