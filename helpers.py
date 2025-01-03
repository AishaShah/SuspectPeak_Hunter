# helpers.py
import pandas as pd

def get_samples(st):
    """
    Retrieve a list of sample IDs from the index of the DataFrame i.e rownames.
    
    Parameters:
    - st (pd.DataFrame): DataFrame representing the samplesheet.
    
    Returns:
    - List of sample IDs.
    """
    return list(st.index)

def get_control(sample_id, st):
    """
    Retrieve the control sample ID for a given sample.
    
    Parameters:
    - sample_id (str): ID of the sample.
    - st (pd.DataFrame): DataFrame representing the samplesheet, 
      with columns 'sample' and 'control'.
    
    Returns:
    - str: Control sample ID for the given sample.
    """
    return st.loc[st['sample'] == sample_id, 'control'].values[0]

def get_type(sample_id, st):
    """
    Retrieve the library type (e.g., Paired End or Single End) for a given sample.
    
    Parameters:
    - sample_id (str): ID of the sample.
    - st (pd.DataFrame): DataFrame representing the samplesheet, with a column 'type'.
    
    Returns:
    - str: Library type for the given sample.
    """
    return st.loc[st['sample'] == sample_id, 'type'].values[0]

def get_Target_Group(sample_id, st):
    """
    Determine the target group for a sample (histone or factor).
    
    Parameters:
    - sample_id (str): ID of the sample.
    - st (pd.DataFrame): DataFrame representing the samplesheet, with a column 'Target_Group'.
    
    Returns:
    - str: 'histone' if target group is active_mark/inactive_mark,
           'factor' if target group is TF,
           or an error message if the value is invalid.
    """
    target_group = st.loc[st['sample'] == sample_id, 'Target_Group'].values[0]
    histone_marks = ["active_mark", "inactive_mark"]
    
    if target_group in histone_marks:
        return "histone"
    elif target_group == "TF":
        return "factor"
    
    return "Please check your metadata file. Target_Group should be 'active_mark', 'inactive_mark', or 'TF'."

def get_sampleIDs_with_IgG(st):
    """
    Retrieve a list of sample IDs that have control samples (IgG).
    
    Parameters:
    - st (pd.DataFrame): DataFrame representing the samplesheet, with a column 'control'.

    Returns:
    - list: List of sample IDs with control samples.
    """
    samples_with_IgG = []
    for sample in st['sample']:
        if pd.notna(st.loc[st['sample'] == sample, 'control'].values[0]):
            samples_with_IgG.append(sample)
    return samples_with_IgG

def IgG_Present(sample, st):
    """
    Check if a control sample exists for a given sample, and retrieve its ID.
    
    Parameters:
    - sample (str): ID of the sample.
    - st (pd.DataFrame): DataFrame representing the samplesheet, with a column 'control'.
    
    Returns:
    - str: Control sample ID if it exists, otherwise 'none'.
    """
    if 'control' in st.columns:
        if pd.notna(st.loc[st['sample'] == sample, 'control'].values[0]):
            return get_control(sample, st)
    return "none"


def get_outfiles(test_samps=None, test_samps_st=None, type="all", array=None):
    """
    return a list of samples PE, SE or both from a given array file of or pandas dataframe
    """
    outfiles_SE = []
    outfiles_PE = []
    if array:
        # read samplesheet and set "sample" name column as row index
        test_samps_st = pd.read_table(array).set_index('sample', drop=False)
        test_samps = get_samples(test_samps_st)  # get sample names
    
    for sample in test_samps:
        sample_type = get_type(sample, test_samps_st)
        
        if sample_type == "PAIRED":
            outfiles_PE.append(sample)
        elif sample_type == "SINGLE":
            outfiles_SE.append(sample)
    
    if type == "SE":
        return list(outfiles_SE)
    elif type == "PE":
        return list(outfiles_PE)
    elif type == "all":
        return list(outfiles_SE + outfiles_PE)
    else:
        raise ValueError("Invalid type parameter. Must be 'SE', 'PE', or 'all'.")


def get_outfiles(test_samps=None, test_samps_st=None, type="all", array=None):
    """
    Retrieve a list of samples based on library type (PE, SE, or both).
    
    Parameters:
    - test_samps (list, optional): List of sample IDs to process. Default is None.
    - test_samps_st (pd.DataFrame, optional): DataFrame representing the samplesheet.
      If not provided, it will be generated from the array file. Default is None.
    - type (str): The type of samples to retrieve. 
      Options are:
      - "SE": Single End samples.
      - "PE": Paired End samples.
      - "all": Both SE and PE samples. Default is "all".
    - array (str, optional): Path to the file containing the samplesheet.
      If provided, this file will be used to generate `test_samps_st`. Default is None.
    
    Returns:
    - list: List of sample IDs matching the specified type (SE, PE, or all).
    
    Raises:
    - ValueError: If the `type` parameter is not one of "SE", "PE", or "all".
    """
    outfiles_SE = []  # List to store Single End sample IDs.
    outfiles_PE = []  # List to store Paired End sample IDs.

    if array:
        # Load the samplesheet from the provided file path and set the "sample" column as the index.
        test_samps_st = pd.read_table(array).set_index('sample', drop=False)
        test_samps = get_samples(test_samps_st)  # Retrieve the list of sample names.

    for sample in test_samps:
        sample_type = get_type(sample, test_samps_st)  # Determine the type of the sample.
        
        if sample_type == "PAIRED":
            outfiles_PE.append(sample)  # Add to Paired End list.
        elif sample_type == "SINGLE":
            outfiles_SE.append(sample)  # Add to Single End list.

    if type == "SE":
        return list(outfiles_SE)
    elif type == "PE":
        return list(outfiles_PE)
    elif type == "all":
        return list(outfiles_SE + outfiles_PE)
    else:
        raise ValueError("Invalid type parameter. Must be 'SE', 'PE', or 'all'.")


def list_array_files(round, path_prefix="", path_inter="", path_suffix="", array_dir=""):
    """
    Generate a list of file paths for samples in specified bootstrap round.

    Parameters:
    - round (int or range): Single bootstrap round number or a range of rounds.
    - path_prefix (str): Prefix path to prepend to the generated file paths. Default is an empty string.
    - path_inter (str): Intermediate path segment to include in the file paths. Default is an empty string.
    - path_suffix (str): Suffix to append to the sample names (e.g., file extensions). Default is an empty string.
    - array_dir (str): Directory containing the bootstrap array files having sample metadata. Default is an empty string.

    Returns:
    - list: A list of file paths generated for the samples in the specified rounds.

    Raises:
    - FileNotFoundError: If the required array file does not exist.

    Example:
    - Input: round=1, path_prefix="data", path_inter="results", path_suffix=".fastq", array_dir="arrays"
    - Output: ["data/round_1/results/sample1.fastq", "data/round_1/results/sample2.fastq", ...]
    """
    import os
    import pandas as pd

    # Ensure `round` is iterable and i.e list format
    rounds = round if isinstance(round, range) else [round]
    #if isinstance(round, range):
    #    rounds = round
    #else:
    #    rounds = [round]

    all_files = []
    for rnd in rounds:
        # Construct the path to the array file
        array = f"{array_dir}/array_{rnd}.tsv"
        if not os.path.exists(array):
            raise FileNotFoundError(f"Required file {array} not found. Ensure bootstrap_samples rule has completed.")
        
        # Read the sample column from the array file
        samples = pd.read_table(array)['sample'].tolist()  # Adjust the column name if needed
        
        # Generate file paths for each sample
        files = [f"{path_prefix}/round_{rnd}/{path_inter}/{sample}{path_suffix}" for sample in samples]
        all_files.extend(files)
    
    return all_files

def check_if_adaptors_provided(st, sample, trim_param):
    """
    Check if forward and reverse adaptors are provided, and generate a Trim Galore command line to use provided adaptores.

    Parameters:
    - st (pd.DataFrame): DataFrame containing the sample metadata.
    - sample (str): Sample ID for which to check adaptors.
    - trim_param (str): Existing Trim Galore parameters to append adaptors to.

    Returns:
    - str: Updated Trim Galore parameters including adaptor sequences if provided.

    Notes:
    - The `Fwd_adaptor` column should contain the forward adaptor sequence.
    - The `Rev_adaptor` column should contain the reverse adaptor sequence (for paired-end libraries).
    - Uses the `get_type()` function to determine if the sample is paired-end.

    Example:
    - Input: trim_param="", sample="sample1"
    - Output: "-a 'FWD_ADAPTOR_SEQ' -a2 'REV_ADAPTOR_SEQ'"
    """
    if pd.notna(st.loc[st['sample'] == sample, 'Fwd_adaptor'].values[0]):
        # Append forward adaptor
        fwd_adaptor = st.loc[st['sample'] == sample, 'Fwd_adaptor'].values[0]
        trim_param = f"{trim_param} -a '{fwd_adaptor}'"
    
    if pd.notna(st.loc[st['sample'] == sample, 'Rev_adaptor'].values[0]):
        # Append reverse adaptor if sample is paired-end
        if get_type(sample, st) == "PAIRED":
            rev_adaptor = st.loc[st['sample'] == sample, 'Rev_adaptor'].values[0]
            trim_param = f"{trim_param} -a2 '{rev_adaptor}'"
    
    return trim_param
