# Configuration Guide for SuspectPeak_Hunter

This guide details the `config.yaml` settings for the **SuspectPeak_Hunter** Snakemake pipeline. Each section outlines the default parameters. [, their purpose]



## General Settings

### `workdir`

The working directory where all files, including Snakemake code, input files, intermediate results, and outputs, will be stored.

**Default**: `/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter`

Set this to a directory that best suits your organizational or storage needs.



## Genome Information

### `genome`

Specifies the genome assembly file for the analysis.

**Default**: `Danio_rerio.GRCz11.dna_sm.primary_assembly`



## Sample Information

### `SL_array`

The `.tsv` file containing metadata for samples used in suspect list generation.

**Default**: `samplesheet.tsv`

Replace this with your custom sample list file.

### `ApplySL_array` [THIS RULE HAS BEEN REMOVED]

The `.tsv` file with metadata for samples where the suspect list will be applied. If `apply_SL` is set to `True`, suspect regions identified in `SL_array` will be excluded from these samples before peak calling.

**Default**: `None`

[THIS SHOULD BE UPDATED --\> section ApplySL in config.yaml having parameters 1) array: samples.tsv and 2) apply_SL:True or False]



## File Settings

These settings determine whether intermediate files are stored or removed.\
Adjust them based on your disk space and the need for downstream analysis.

-   **`keep_trimmed_reads_fastq`**: `True`
-   **`keep_Raw_BAMs`**: `False`
-   **`keep_BAMs_wo_MT`**: `True` BAM files after removing Mitochondrial Genome
-   **`keep_filtered_BAMs`**: `True` [BAM files after rmeoving mitochondrial genome and multimapped reads [HOW ARE MM READS REMOVED \<-- CONFIRM \<-- Removing reads having "XS" tag that stores information of secondary alignment \<-- maybe we can allow user to choose if they want to remove all MM reads or keep primary alignments of MM reads and just remove secondary]]{style="color:red;"}
-   **`keep_genomecov_bedgraph`**: `True`
-   **`keep_Peaks`**: `True`



## Trimming and Mapping

### Trimming Settings

#### `skip_trimming`

Set to `True` to skip trimming if reads are already processed.

**Default**: `False`

#### Trimgalore Parameters

Adjust trimming stringency, quality thresholds, and threads as needed.

-   **Paired-end**:

    -   **Additional Params**: `--quality 20 --length 20 --stringency 1 -e 0.1`
    -   **Threads**: `24`

-   **Single-end**:

    -   **Additional Params**: `--illumina -q 20 --length 20 --stringency 5 -e 0.5`
    -   **Threads**: `24`

### Mapping Settings

#### Bowtie2 Parameters

Modify these settings to fine-tune alignment sensitivity and specificity.

-   **Single-end**:

    -   **Additional Params**: `--very-sensitive-local --no-unal`
    -   **Threads**: `24`

-   **Paired-end**:

    -   **Additional Params**: `--very-sensitive-local --no-mixed --no-unal --dovetail -X 1000`
    -   **Threads**: `24`



## Downsampling

### `DownSample`

Adjust read counts and threads to optimize for dataset size and computational constraints. We suggest downsampling so that samples are more comparable and less computation resources will be used.

-   **Reads**: `10000000`
-   **Reads_PE**: `20000000`
-   **Reads_SE**: `10000000`
-   **Threads**: `8`



## Peak Calling

### `peak_calling`

Refine the sensitivity and specificity of peak detection using the following parameters:

-   **Mode**: `stringent`
-   **Threshold**: `0.1`
-   **Normalization**: `non`
-   **Filter by Depth**: `True`
-   **Min Depth**: `5`



# Parameters for bootstrapping

## Suspect List Generation

### `bootstrap_samples`

Defines the parameters for bootstrapping, which helps improve the robustness of suspect region identification by iterating through random subsets of the data:

-   **Rounds**: The number of bootstrap iterations.

-   **Samples Per Group**: The number of samples included in each sample group. Groups are defined by column named "Target_Groups" in input samplesheet.tsv file i.e SL_array. For example in our case there were three groups i.e TF, Active marks, Inactive marks. We try to take equal number of sample from each group in each bootstrap iteration.

-   **Seed**: The random seed used for reproducibility of the bootstrap sampling process.

**Default Values**:

-   **Bootstrap Rounds**: `10`

-   **Samples Per Group**: `3`

-   **Seed**: `1`

### `generate_suspectList`

This section defines the parameters for identifying suspect regions across samples in each bootstrap iteration. The suspect list is generated based on the following criteria:

-   **Number of Samples**: Specifies the number of samples randomly selected during each bootstrap iteration.

-   **Minimum Length**: The minimum length (in base pairs) of a peak region for it to be considered for inclusion in the suspect list.

-   **Percentage Threshold**: The minimum percentage of samples in which a peak region must appear to qualify as a suspect region.

**Default Values**:

-   Number of Samples: `9`
-   Minimum Length: `10`
-   Percentage Threshold: `80`

Adjust these parameters based on the characteristics of your dataset to optimize suspect region identification.



## MACS2 Peak Calling

### `macs2_callpeak`

Parameters for calling peaks using MACS2.

-   **Mappable Genome Size**: `1368780147`

Update this value to match the genome size of your organism for more accurate peak calling.\

------------------------------------------------------------------------

# Example of full config.yaml file:

``` yaml
workdir: "/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter"
# Genome information
genome: "Danio_rerio.GRCz11.dna_sm.primary_assembly"

# Sample information
SL_array: "array_rep1.wo_IgG.wo_NegCtrls.tsv"
ApplySL_array: "PRJNA738523.samples.array"

# File settings
keep_trimmed_reads_fastq: "True"
keep_Raw_BAMs: "False"
keep_BAMs_wo_MT: "True"
keep_filtered_BAMs: "True"
keep_genomecov_bedgraph: "True"
keep_Peaks: "True"

# Trimming and mapping settings
skip_trimming: "False"  # Set to "True" to skip trimming
skip_mapping: "False"   # Set to "True" to skip mapping

# Input file type provided:
# BAMs_provided: "False"
# Trimmed_fastq_provided: "False"
# Fastq_provided: "True"

# Parameters for Trimgalore
trim_galore_pe:
  additional_params: "--quality 20 --length 20 --stringency 1 -e 0.1"
  threads: 24
trim_galore_se:
  additional_params: "--illumina -q 20 --length 20 --stringency 5 -e 0.5"
  threads: 24

# Mapping parameters for bowtie2
map_bowtie2_se:
  additional_params: "--very-sensitive-local --no-unal"
  threads: 24
map_bowtie2_pe:
  additional_params: "--very-sensitive-local --no-mixed --no-unal --dovetail -X 1000"
  threads: 24

bootstrap_samples:
  rounds: 10
  samples_per_group: 3
  seed: 1

DownSample:
  Reads: "10000000" #"3000000" 
  Reads_PE: "20000000" #"6000000"
  Reads_SE: "10000000" #"3000000"
  threads: 8
GenMap:
  MinMap: "0.3"

# Peak calling settings
peak_calling:
  mode: "stringent"
  threshold: "0.1"  # "0.1" "0.2" # "0.5" "0.9" "1" "0.001"
  normalization: "non"
  filter_by_depth: True
  min_depth: 5

# Suspect list generation settings
generate_suspectList:
  num_of_samples: "9"
  minimum_length: "10"
  percentage_threshold: "80"


# Peak calling using MACS2
macs2_callpeak:
  mappable_genome_size: 1368780147
```

------------------------------------------------------------------------

## Additional Notes

-   Default parameters are chosen to optimize performance for standard datasets.

-   Modify parameters in `config.yaml` to suit your specific experimental setup and computational resources.

-   Future updates will include more configuration options and automation for enhanced flexibility.
