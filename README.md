# SuspectPeak_Hunter

## Vastenhouw Lab - Aisha Shah

### **Overview**

SuspectPeak_Hunter is a Snakemake pipeline designed to enhance CUT&RUN data analysis by identifying "suspect regions" in the genome. These regions are consistently present across multiple samples and often appear in the top-ranked peaks, potentially overshadowing more biologically relevant, target-specific peaks. By providing a list of these recurring regions, SuspectPeak_Hunter enables users to exclude them before peak calling, ensuring that their analyses focus on peaks more likely to be relevant to the histone marks or transcription factors under study. This approach helps improve the accuracy and interpretability of CUT&RUN results, particularly for workflows relying on top-ranked peaks. This pipeline automates various steps from quality control, read trimming, genome mapping, to peak calling and finally suspect list generation. (THIS WAS IMPLEMENTED IN PREVIOUS VERSION) It can also be used to apply generated suspect list and call peaks after removing suspect-list regions on a given set of samples.

### **Main Steps:**

-   Quality control of raw reads using FastQC
-   Read trimming using Trim Galore
-   Genome indexing and mapping using Bowtie2
-   Peak calling with SEACR
-   Suspect list generation by identifying regions present across samples in different groups (i.e TF, Active marks, inactive marks)
-   Call Peaks after applying suspect list (THIS WAS REMOVED FROM CURRENT VERSION)

### **Installation**

**Clone the repository**: 

```sh
git clone https://github.com/yourusername/SuspectPeak_Hunter.git     
cd SuspectPeak_Hunter
```

**Option A: Installing in Conda environment**

1.  Create and activate a conda environment: 

      ```sh     
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

1. Download the samples.
2. Prepare the samplesheet.
3. Set up the configuration file.


**Downloading Samples** 

Samples can be downloaded either manually by the user or using the provided script: `scripts/download_samples.py`. This script requires an input file containing SRR sample IDs (one ID per line).

Example Input File:
examples/ids.txt

```yaml
SRR328768  
SRR865977  
ERR099467  
```
    
Command to Run the Script:
```sh
python scripts/download_samples.py examples/ids.txt   
```
Replace examples/ids.txt with the path to your input file containing sample IDs.


**Prepare the sample sheet** **(`samplesheet.tsv`)**

A tab-separated file containing the same columns as the example. The columns of the sample file are:

-   **projectID**

-   **ctrl:** is this a control sample or target enriched

-   **shtname:** optional \<\-- not required \--\> remove

-   **sample:** Sample name

-   **control:** name of control sample if exists

-   **rep:** replicate id of the sample. this is required for during bootstraping we not allow replicates of same sample in same bootstrap round

-   **type:** PAIRED or SINGLE. required to choose mapping configuration and genomecove file \
    configuration path: path to fastq file

-   **R1:** name of read1/fwd pairs fastq

-   **R2:** name of read2/rev pairs fastq (if type==PAIRED)


| projectID | ctrl            | shtname | sample | control | rep | type   | path           | R1         | R2         |
|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|
| PRJXX     | target_enriched | sample1 | SRRXX  | SRRXX   | 1   | PAIRED | /path/to/fastq | SRRXX_1.fq | SRRXX_2.fq |
| PRJXX     | target_enriched | sample1 | SRRXX  | SRRXX   | 2   | PAIRED | /path/to/fastq | SRRXX_1.fq | SRRXX_2.fq |
| PRJXX     | neg_ctrl        | sample2 | SRRXX  | SRRXX   | 1   | SINGLE | /path/to/fastq | SRRXX.fq   | NA         |
| PRJXX     | target_enriched | sample1 | SRRXX  | NA      | 1   | PAIRED | /path/to/fastq | SRRXX_1.fq | SRRXX_2.fq |



**Prepare the configuration file (`config.yaml`)**:

(More input options to be added in future)

```yaml
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

### **Run the pipeline**:

-   Run locally: `sh       snakemake --cores <number_of_cores>`

-   Run on HPC cluster: `sh       snakemake XX`

-   Run using docker: `sh       snakemake XX`

### Pipeline Steps

#### Quality Control

Perform quality control of raw reads using FastQC.

#### Read Trimming

Trim reads using Trim Galore to remove adapters and low-quality bases.

#### Genome Indexing

Index the genome using Bowtie2 for efficient mapping.

#### Mapping

Map the trimmed reads to the genome using Bowtie2.

#### Peak Calling

Identify peaks in the mapped data using SEACR.

#### Generating Suspect Lists

Generate a list of suspect regions in the genome based on the peak calling results.

### Contribution

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
