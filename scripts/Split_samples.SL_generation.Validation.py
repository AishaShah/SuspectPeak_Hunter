
import pandas as pd
import numpy as np
import argparse
from collections import defaultdict
import os


def load_metadata(file_path):
    return pd.read_csv(file_path, sep='\t')

def balance_split(metadata, test_size=0.3, seed=None, replicate_id=None):
    if seed is not None:
        np.random.seed(seed)

    # Filter by replicate ID if provided
    if replicate_id:
        metadata = metadata[metadata['rep'] == replicate_id]

    df_for_target_rarity = metadata[['shtname', 'Target', 'Stage']].drop_duplicates()

    # Count occurrences of each target
    target_counts = df_for_target_rarity['Target'].value_counts()

    # Identify rare targets (those that occur only once)
    rare_targets = target_counts[target_counts == 1].index

    # Exclude rare targets from test set consideration
    eligible_targets = metadata[~metadata['Target'].isin(rare_targets)]

    # Get unique combinations of Target and Stage for eligible targets
    eligible_combinations = eligible_targets[['Target', 'Stage','Target_Group']].drop_duplicates()

    # Ensure at least one sample from each target group in the test set
    target_groups = metadata['Target_Group'].unique()  # Get all unique target groups from original metadata
    
    test_targets = pd.DataFrame()

    for group in target_groups:
        group_combinations = eligible_combinations[eligible_combinations['Target_Group'] == group]
        if group_combinations.empty:
            # Print a warning if the target group is not present in eligible combinations
            print(f"Warning: Target group '{group}' is not included in the test set as it's not present in the eligible dataset.")
        else:
            # Ensure at least one combination per group in the test set
            selected = group_combinations.sample(n=1, random_state=seed)
            test_targets = pd.concat([test_targets, selected])

    # Sample the remaining combinations based on the test size
    remaining_combinations = eligible_combinations[~eligible_combinations.index.isin(test_targets.index)]
    
    additional_test_targets = remaining_combinations.sample(
        frac=(test_size - (len(test_targets) / len(metadata[['Target', 'Stage','Target_Group']].drop_duplicates()))),
        random_state=seed
    )
    #metadata_for_additional_test_targets=metadata[metadata[['Target', 'Stage']].apply(tuple, axis=1).isin(additional_test_targets.apply(tuple, axis=1))]
    test_targets = pd.concat([test_targets, additional_test_targets])

    # Include all rows with the selected targets and stages in the test set
    test_metadata = metadata[metadata[['Target', 'Stage', 'Target_Group']].apply(tuple, axis=1).isin(test_targets.apply(tuple, axis=1))]

    # The remaining data will be the training set
    train_metadata = metadata[~metadata.index.isin(test_metadata.index)]

    return train_metadata, test_metadata



def parse_args():
    parser = argparse.ArgumentParser(description='Split sample metadata into SL_generation_samples and validation_samples.')
    parser.add_argument('metadata_file', type=str, help='Path to the metadata file.')
    parser.add_argument('--validation_set_size', type=float, default=0.3, help='Proportion of the datasets to include in the validation set.')
    parser.add_argument('--replicate_id', type=str, default=None, help='Filter metadata by replicate ID.') #### NOT NEEDED
    parser.add_argument('--seed', type=int, default=None, help='Seed for random number generation.')
    parser.add_argument('--output_path', type=str, default=None, help='Directory for results.')

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    metadata = load_metadata(args.metadata_file)
    SL_generation_samples, validation_samples  = balance_split(
        metadata,
        test_size=args.validation_set_size,
        seed=args.seed,
        replicate_id=args.replicate_id
    )
    
    output_dir = args.output_path or os.getcwd()
    os.makedirs(output_dir, exist_ok=True)
    SL_generation_samples_output = os.path.join(output_dir, 'SL_generation_samples.tsv')
    validation_samples_output = os.path.join(output_dir, 'validation_samples.tsv')

    SL_generation_samples.to_csv(SL_generation_samples_output, sep='\t', index=False,na_rep='NA')
    validation_samples.to_csv(validation_samples_output, sep='\t', index=False,na_rep='NA')