# SuspectPeak_Hunter

## Vastenhouw Lab - Aisha Shah

### Overview
SuspectPeak_Hunter is a Snakemake pipeline designed to generate a suspect list of problematic regions in the genome for peak calling analysis. This pipeline automates various steps from quality control, read trimming, genome mapping, to peak calling and finally suspect list generation. This pipline can also be use to apply generated suspect list and call peaks after removing suspect-list regions on a given set of samples.

### Features
- Quality control of raw reads using FastQC
- Read trimming using Trim Galore
- Genome indexing and mapping using Bowtie2
- Peak calling with SEACR
- Suspect list generation by identifying problematic regions in the genome
- Call Peaks after applying suspect list

### Installation
1. **Clone the repository**:
    ```sh
    git clone https://github.com/yourusername/SuspectPeak_Hunter.git
    cd SuspectPeak_Hunter
    ```
2. **Install required software**:
    - Snakemake
    - FastQC
    - Trim Galore
    - Bowtie2
    - Samtools
    - Bedtools
    - SEACR
    - R

3. **Create and activate a conda environment** (optional but recommended):
   (** TO BE UPDATED **)
    ```sh
    conda create -n suspectpeak_hunter python=3.8
    conda activate suspectpeak_hunter
    pip install -r requirements.txt
    ```

### Usage
1. **Prepare the configuration file (`config.yaml`)**:
    (** More input options to be added in future **)
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

2. **Prepare the sample sheet (`samplesheet.tsv`)**:
(TO BE UPDATED)  

| projectID | ctrl            | shtname      | sample | control | rep | type   | path            | R1         | R2         |
|-----------|-----------------|--------------|--------|---------|-----|--------|-----------------|------------|------------|
| PRJXX     | target_enriched | sample1      | SRRXX  | SRRXX   | 1   | PAIRED | /path/to/fastq  | SRRXX_1.fq | SRRXX_2.fq |
| PRJXX     | target_enriched | sample1      | SRRXX  | SRRXX   | 2   | PAIRED | /path/to/fastq  | SRRXX_1.fq | SRRXX_2.fq |
| PRJXX     | neg_ctrl        | sample2      | SRRXX  | SRRXX   | 1   | SINGLE | /path/to/fastq  | SRRXX.fq   |  NA        |
| PRJXX     | target_enriched | sample1      | SRRXX  | NA      | 1   | PAIRED | /path/to/fastq  | SRRXX_1.fq | SRRXX_2.fq |


3. **Run the pipeline**:

  -  Run locally:
      ```sh
      snakemake --cores <number_of_cores>
      ```
  -  Run on HPC cluster:
      ```sh
      snakemake XX
      ```
    
  -  Run using docker:
      ```sh
      snakemake XX
      ```
  

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
1. Fork the repository
2. Create a new branch (`git checkout -b feature-branch`)
3. Commit your changes (`git commit -am 'Add new feature'`)
4. Push to the branch (`git push origin feature-branch`)
5. Open a pull request

### License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

### Acknowledgements
This pipeline was developed by Aisha Shah at the Vastenhouw Lab under the supervision of Prof. Nadine Vastenhouw (PI) and Dora Grabavac (PhD).

---

For any issues or questions, please contact Aisha Shah at [aisha.shah@alumni.esci.upf.edu](mailto:aisha.shah@alumni.esci.upf.edu).
