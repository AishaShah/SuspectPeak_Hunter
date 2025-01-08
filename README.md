# SuspectPeak_Hunter

## Vastenhouw Lab - Aisha Shah

### Overview

SuspectPeak_Hunter is a Snakemake pipeline designed to enhance CUT&RUN data analysis by identifying "suspect regions" in the genome. These regions are consistently present across multiple samples and often appear in the top-ranked peaks, potentially overshadowing more biologically relevant, target-specific peaks. By providing a list of these recurring regions, SuspectPeak_Hunter enables users to exclude them before peak calling, ensuring that their analyses focus on peaks more likely to be relevant to the histone marks or transcription factors under study. This approach helps improve the accuracy and interpretability of CUT&RUN results, particularly for workflows relying on top-ranked peaks. This pipeline automates various steps from quality control, read trimming, genome mapping, to peak calling and finally suspect list generation. (THIS WAS IMPLEMENTED IN PREVIOUS VERSION) It can also be used to apply generated suspect list and call peaks after removing suspect-list regions on a given set of samples.

### Main Steps:

-   Quality control of raw reads using FastQC
-   Read trimming using Trim Galore to remove adapters and low-quality bases
-   Genome indexing and mapping using Bowtie2
-   Peak calling with SEACR
-   Suspect list generation by identifying regions present across samples in different groups (i.e TF, Active marks, inactive marks)
-   Call Peaks after applying suspect list (THIS WAS REMOVED FROM CURRENT VERSION)

### Installation

**Clone the repository**:

``` sh
git clone https://github.com/yourusername/SuspectPeak_Hunter.git     
cd SuspectPeak_Hunter
```

**Option A: Installing in Conda environment**

1.  Create and activate a conda environment:

    ``` sh
    suspeak_hunter_env="/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/conda_env/suspeak_hunter_env_mamba"   
    conda_pkgs_dir="/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/conda_pkgs_dirs/"     
    mamba env create --prefix $suspeak_hunter_env -f trimgalore_env.yml     
    export CONDA_PKGS_DIRS=$conda_pkgs_dir     
    mamba install -c conda-forge pandas     
    conda activate $suspeak_hunter_env
    ```

2.  Install required software from official websites:

    -   [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
    -   [FastQC](https://anaconda.org/bioconda/fastqc)
    -   [Trim Galore](https://github.com/FelixKrueger/TrimGalore)
    -   [Bowtie2](https://github.com/BenLangmead/bowtie2)
    -   [Samtools](https://www.htslib.org/download/)
    -   [Bedtools](https://www.htslib.org/download/)
    -   [SEACR](https://github.com/FredHutch/SEACR)
    -   [R](https://rstudio-education.github.io/hopr/starting.html)

**Option B: Installing Docker**

[DOCKER SECTION TO BE UPDATED]

### Usage

Before running the pipeline, you need to:

1.  Download the samples.
2.  Prepare the samplesheet.
3.  Set up the configuration file.

[**Downloading Samples**]{.underline}

Samples can be downloaded either manually by the user or using the provided script: `scripts/download_samples.py`. This script requires an input file containing SRR sample IDs (one ID per line).

Example Input File: examples/ids.txt

``` yaml
SRR328768  
SRR865977  
ERR099467  
```

Command to Run the Script:

``` sh
python scripts/download_samples.py examples/ids.txt   
```

Replace examples/ids.txt with the path to your input file containing sample IDs.

[**Prepare the sample sheet(`samplesheet.tsv`)**]{.underline}

A tab-separated file containing the same columns as the example. The columns of the sample file are:

-   **projectID**

-   **ctrl:** is this a control sample or target enriched

-   **shtname:** optional \<\-- not required \--\> remove

-   **sample:** Sample name

-   **control:** name of control sample if exists

-   **rep:** replicate id of the sample. this is required for during bootstraping we not allow replicates of same sample in same bootstrap round

-   **type:** PAIRED or SINGLE. required to choose mapping configuration and genomecove file\
    configuration path: path to fastq file

-   **R1:** name of read1/fwd pairs fastq

-   **R2:** name of read2/rev pairs fastq (if type==PAIRED)

| projectID | Group         | ctrl            | shtname | sample | control | rep | type   | path           | R1         | R2         |
|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
| PRJXX     | TF            | target_enriched | sample1 | SRRXX  | SRRXX   | 1   | PAIRED | /path/to/fastq | SRRXX_1.fq | SRRXX_2.fq |
| PRJXX     | Active Mark   | target_enriched | sample1 | SRRXX  | SRRXX   | 2   | PAIRED | /path/to/fastq | SRRXX_1.fq | SRRXX_2.fq |
| PRJXX     | Active Mark   | neg_ctrl        | sample2 | SRRXX  | SRRXX   | 1   | SINGLE | /path/to/fastq | SRRXX.fq   | NA         |
| PRJXX     | Inactive Mark | target_enriched | sample1 | SRRXX  | NA      | 1   | PAIRED | /path/to/fastq | SRRXX_1.fq | SRRXX_2.fq |

To ensure robust analysis, we recommend dividing your datasets into two groups:

1.  **Suspect List Generation**: Use approximately **70%** of your samples to generate the suspect list through bootstrapping.

2.  **Validation**: Use the remaining **30%** of your samples to validate the suspect list by assessing whether the identified suspect regions are present in these samples.

The above mentioned `samplesheet.tsv` file, referenced in the `config.yaml`, contains metadata for all samples used in the pipeline. By default, all samples in the file are utilized for suspect list generation. To split your samples, use the provided script `scripts/split_samplesheet_data.py` which you can run using the following commandline:

``` python
python scripts/split_samplesheet_data.py --input samplesheet.tsv --output_dir split_samples
```

The script will generate two output files:

-   `samplesheet_train.tsv`: For suspect list generation.

-   `samplesheet_validate.tsv`: For validation.

The script uses the **`Group`** column in `samplesheet.tsv` to stratify samples. It ensures **30% of samples** are selected proportionally from each group for validation. Remaining samples can be used for suspect list generation.

[**Prepare the configuration file (`config.yaml`)**:]{.underline}

(More input options to be added in future)

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

# Suspect list generation settings
generate_suspectList:
  num_of_samples: "9"
  minimum_length: "10"
  percentage_threshold: "80"
```

For detailed information of config.yaml file please check: [wiki/CONFIGURATION.md](https://github.com/AishaShah/SuspectPeak_Hunter/blob/main/wiki/CONFIGURATION.md)

### Running the pipeline

-   Run locally: `snakemake --cores <number_of_cores> Initialize_SuspectPeak_Hunter`

-   Run on HPC cluster: `snakemake --cores <number_of_cores> Initialize_SuspectPeak_Hunter --cluster cluster.yaml "{cluster.nodes} {cluster.time} ..."`

    The file cluster.yaml contain configuration for running on cluster. for each rule we can define nodes, time etc. if not specified for a rule then default parameters given by **default** would be used.

    ``` yaml
    rule A:
      nodes: 10
      time:00:10:00
    ```

-   Run using docker: `snakemake XX`

### To Contribute

1.  Fork the repository
2.  Create a new branch (`git checkout -b feature-branch`)
3.  Commit your changes (`git commit -am 'Add new feature'`)
4.  Push to the branch (`git push origin feature-branch`)
5.  Open a pull request

### License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

### Acknowledgements

This pipeline was developed by Aisha Shah at the Vastenhouw Lab under the supervision of Prof. Nadine Vastenhouw (PI) and Dora Grabavac (PhD).

------------------------------------------------------------------------

For any issues or questions, please contact Aisha Shah at [aisha.shah\@alumni.esci.upf.edu](mailto:aisha.shah@alumni.esci.upf.edu).
