# SuspectPeak_Hunter

## Vastenhouw Lab - Aisha Shah

### Overview
SuspectPeak_Hunter is a Snakemake pipeline designed to generate a suspect list of problematic regions in the genome for peak calling analysis. This pipeline automates various steps from quality control, read trimming, genome mapping, to peak calling and suspect list generation.

### Features
- Quality control of raw reads using FastQC
- Read trimming using Trim Galore
- Genome indexing and mapping using Bowtie2
- Peak calling with SEACR
- Suspect list generation by identifying problematic regions in the genome

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
    genome: "Danio_rerio.GRCz11.dna_sm.primary_assembly"
    peak_calling:
      mode: "stringent"
      threshold: "0.001"
      normalization: "non"
    generate_suspectList:
      percentage_threshold: "60"
      num_of_samples: "32"
      minimum_length: "10"
    ```

2. **Prepare the sample sheet (`samplesheet.tsv`)**:
(TO BE UPDATED)
    ```tsv
    sample    control    type
    sample1   control1   PAIRED
    sample2   control2   SINGLE
    ```

3. **Run the pipeline**:
    ```sh
    snakemake --cores <number_of_cores>
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
Analysis
