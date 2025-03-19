# File path: generate_configs.py
import os
import itertools
import yaml
import argparse
import sys

def read_base_config(file_path):
    """Read the base config.yaml file."""
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)

def collect_available_keys(config, prefix=""):
    """Recursively collect all keys from the YAML config."""
    keys = []
    for key, value in config.items():
        full_key = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            keys.extend(collect_available_keys(value, full_key))
        else:
            keys.append(full_key)
    return keys

def check_key_exists(config, key):
    """Check if a dot-separated key exists in the nested dictionary."""
    keys = key.split('.')
    temp = config
    for k in keys:
        if k not in temp:
            raise KeyError(f"Error: Parameter '{key}' does not exist in the base configuration.")
        temp = temp[k]

def generate_configs(base_config, params, output_dir):
    """Generate combinations of config files based on the parameters."""
    # Validate all parameters
    for param_key in params.keys():
        check_key_exists(base_config, param_key)
    
    # Generate all combinations of parameters
    param_keys = list(params.keys())
    param_values = list(params.values())
    combinations = list(itertools.product(*param_values))
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate a new config file for each combination
    for idx, combo in enumerate(combinations, start=1):
        new_config = base_config.copy()
        
        # Update config with current combination
        for key, value in zip(param_keys, combo):
            keys = key.split('.')
            temp = new_config
            for k in keys[:-1]:
                temp = temp[k]
            temp[keys[-1]] = value
        


        # Update the results_dir parameter dynamically
        results_dir = f"results/config_{idx}"
        new_config['results_dir'] = results_dir  # Ensure results_dir exists in the base config
        
        # Write the updated config to a new file
        output_file = os.path.join(output_dir, f"config_{idx}.yaml")
        with open(output_file, 'w') as file:
            yaml.dump(new_config, file)
        print(f"Generated: {output_file} with results_dir={results_dir}")

def custom_help(config_path=None):
    """Display detailed help with available parameters if --help is passed."""
    help_message = """
    Generate multiple config files by varying parameters.

    Usage:
    python generate_configs.py config.yaml --params key=value1,value2,... [--output-dir OUTPUT_DIR]

    Options:
    config              Path to the base config.yaml file.
    --params            Parameters to modify, provided as key=value or key=value1,value2,...\n
                        Use dot notation for nested keys, e.g., 'DownSample.Reads=100,200'.
    --output-dir        Directory to store generated config files (default: "generated_configs").

    For detailed parameter options, pass config file with --help.
    """
    print(help_message)
    if config_path:
        print("\nAvailable Parameters in config.yaml:")
        base_config = read_base_config(config_path)
        available_keys = collect_available_keys(base_config)
        for key in available_keys:
            print(f"  - {key}")
    sys.exit(0)

def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("config", nargs="?", help="Path to the base config.yaml file.")
    parser.add_argument("--output-dir", default="generated_configs", help="Directory to store generated config files.")
    parser.add_argument("--params", nargs='+', 
                        help=("Parameters to modify, provided as key=value or key=value1,value2,...\n"
                              "Use dot notation for nested keys, e.g., 'DownSample.Reads=100,200'."))
    parser.add_argument("-h", action="store_true", help="Show this help message and exit.")
    parser.add_argument("--help", action="store_true", help="Show detailed help with available parameters.")
    
    args = parser.parse_args()
    
    # Handle -h for basic help
    if args.h:
        parser.print_help()
        sys.exit(0)
    
    # Handle --help for detailed help
    if args.help:
        config_path = args.config if args.config else None
        custom_help(config_path)
    
    # Validate config file
    if not args.config:
        print("Error: Config file is required. Use -h or --help for usage.")
        sys.exit(1)
    
    # Read base config
    base_config = read_base_config(args.config)
    
    # Parse parameters
    params = {}
    if args.params:
        for param in args.params:
            key, values = param.split('=')
            params[key] = values.split(',')
    
    # Generate config files
    try:
        generate_configs(base_config, params, args.output_dir)
    except KeyError as e:
        print(e)
        sys.exit(1)

if __name__ == "__main__":
    main()
