# **SuspectPeak_Hunter Wiki**

This repository contains essential files to guide you in configuring the **SuspectPeak_Hunter** pipeline. The two key documents in this directory are:

-   **CONFIGURATION.md**: Explains how to set up and customize the `config.yaml` file.

-   **SAMPLESHEET.md**: Details the structure and content required for the `samplesheet.tsv` file.

------------------------------------------------------------------------

## Configuring the Pipeline Parameters

The `config.yaml` file is central to the pipeline's functionality. It contains all the parameters needed to control each step of the analysis, such as trimming, mapping, peak calling, and suspect list generation.

### How to Use `config.yaml`

1.  Modify parameters as needed for your specific experiment. You can modify parameters for each step of pipline using this file.

2.  Save the file and ensure it is in the pipeline's working directory.

### Key Sections of `config.yaml`

The file includes configurations for:

-   **General Settings**: Working directory and genome assembly.

-   **Sample Information**: Metadata file for input samples for suspect list generation.

-   **Intermediate Files**: Options to save or discard intermediate files.

-   **Trimming and Mapping**: Parameters for quality control and read alignment.

-   **Peak Calling**: Thresholds and modes for peak detection.

-   **Suspect List Generation**: Settings for identifying problematic genomic regions.

For detailed instructions, refer to **CONFIGURATION.md**.

------------------------------------------------------------------------

## Preparing the Sample Sheet

The `samplesheet.tsv` file provides metadata for all the input samples used in the pipeline. This metadata is essential for ensuring the pipeline processes the correct files and performs the required analyses according to sample type.

### Structure of `samplesheet.tsv`

The file must be in a tab-delimited format with atleast the following columns:

1.  **Sample ID**: A unique identifier for each sample.

2.  **Condition**: A label describing the experimental condition (e.g., treatment, control).

3.  **File Path**: The path to the input file(s) for each sample i.e fastq .

4.  **Read Type**: Specifies if the sample contains paired-end (`PE`) or single-end (`SE`) reads.

### Example `samplesheet.tsv`

``` python
Sample_ID   Condition    File_Path                           Read_Type
Sample1     Control      /data/samples/sample1_R1.fastq.gz  PE
Sample2     Treated      /data/samples/sample2_R1.fastq.gz  SE
Sample3     Control      /data/samples/sample3_R1.fastq.gz  PE
```

### Guidelines

-   Ensure all paths are correct and accessible from the working directory.

-   Update the `config.yaml` file to reference your `samplesheet.tsv`.

For more details on formatting and required fields and how these fields are used by the pipline, refer to **SAMPLESHEET.md**.

Both `config.yaml` and `samplesheet.tsv` are required to run the pipeline.
