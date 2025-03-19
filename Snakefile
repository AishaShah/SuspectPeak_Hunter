
######################################################################################################################
#                                                                                                                    #
#                                                SuspectPeakHunter                                                   #
#                                                                                                                    #
#                       A Snakemake pipeline for generating suspect lists for CUT&RUN data                           #
#                                                                                                                    #
######################################################################################################################
                                                                                                                
# path to the config file and in/  and in/out directories 
configfile: 'config.yaml'
workdir: config["workdir"]
resultdir: config["results_dir"]

# load required python packages
import pandas as pd

#################################################
#          Helper Functions for Metadata        #
#################################################

from helpers import *

#*************#
# INPUT FILES #
#*************#
genome = config["genome"] # Load genome


#*****#
# (1) #  Input files for generating suspect list
#*****#

# samplesheet for all samples
all_samples=config["all_samples"] 
# read samplesheet and set "sample" name column as row index
all_samples_df = pd.read_table(all_samples).set_index('sample', drop=False) 


# Load samplesheet for samples needed for suspectlist generation
SL_generation_samples=config["SL_generation_samples"] 
# read samplesheet and set "sample" name column as row index
SL_generation_samples_df = pd.read_table(SL_generation_samples).set_index('sample', drop=False) 
samps = get_samples(SL_generation_samples_df)   # get sample names
seacr_mode=config["peak_calling"]["mode"]
seacr_threshold=config["peak_calling"]["threshold"]


# the following is to make sure that wildcard for sample_name does not contain "." 
# to avoid snakemake considering sample_name=noSL.Sample1 instead of sample_name=Sample1
wildcard_constraints:
    sample="[^.]+"


#*****#
# (2) #  Input files for applying suspect list
#*****#

# single end and paired-end sample names
SE=get_outfiles(samps, SL_generation_samples_df, type="SE")
PE=get_outfiles(samps, SL_generation_samples_df, type="PE")


# The following script generates bootstrap arrays for each round
# The script generates a tsv file for each round containing the samples to be used for that round
# if arrays are already generated, the script will not generate them again
import subprocess
subprocess.run(["python3", "scripts/generate_bootstrap_arrays.py", 
                SL_generation_samples, 
                str(config["bootstrap_samples"]["rounds"]), 
                str(config["bootstrap_samples"]["samples_per_group"]), 
                "04.bootstrapping_arrays", 
                "--seed", str(config["bootstrap_samples"]["seed"])])



#################################################
## Generate SuspectList                        ##
#################################################
# The pipeline generates suspect lists for CUT&RUN data
# The pipeline consists of the following steps:
# (1) Prepare Genome for Mapping
# (2) FastQC and trimming
# (3) Mapping
# (4) Remove mitochondrial chromosome
# (5) BootStrapping
#   (5.1) Downsample samples for each bootstrap round
#   (5.2) Convert BAM files to BED files
#   (5.3) GenomeCov
#   (5.4) Filter GenomeCov by min_depth filter (Optional)
#   (5.5) PeakCalling
#   (5.6) Overlap peaks
#   (5.7) Generate SuspectList
# (6) Overlap SuspectListes from all bootstrap rounds to generate all_rounds.SL.bed
# (7) Validate SuspectList
#   (7.1) Call peaks without removing SL regions
#   (7.1) Remove SuspectList regions from BAMTOBED files
#   (7.2) genomecov
#   (7.3) Call peaks
#   (7.4) Overlap Validation sample peaks with SuspectList



# rule Initialize_SuspectPeak_Hunter runs all the other rules required for generating suspectlist
# The rule generates suspect lists for each bootstrap round and for all rounds
rule Initialize_SuspectPeak_Hunter:
    input:
        # overlap peaks per round
        expand(
            "{result_dir}/05.bootstrapping/05.Overlap_Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.mergepeaks/round_{round}.rawbam.mergedPeaks.binary.tab", 
            result_dir=config["results_dir"],
            num_reads=config["DownSample"]["Reads"],
            round=range(1, config["bootstrap_samples"]["rounds"] + 1),
            #round=range(1, 2 + 1),
            seacr_threshold=config["peak_calling"]["threshold"],
            min_depth=config["peak_calling"]["min_depth"]
         ) if config["peak_calling"]["filter_by_depth"] else  expand(
            "{result_dir}/05.bootstrapping/05.Overlap_Peaks.DS_{num_reads}/SEACR_{seacr_threshold}.mergepeaks/round_{round}.rawbam.mergedPeaks.binary.tab", 
            num_reads=config["DownSample"]["Reads"],
            round=range(1, config["bootstrap_samples"]["rounds"] + 1),
            #round=range(1, 2 + 1),
            seacr_threshold=config["peak_calling"]["threshold"]
            ),
        
        # generate suspect list per round
        expand("{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/round_{round}.SL.bed", 
            result_dir=config["results_dir"],
            num_reads=config["DownSample"]["Reads"],
            round=range(1, config["bootstrap_samples"]["rounds"] + 1),
            #round=range(1, 2 + 1),
            seacr_threshold=config["peak_calling"]["threshold"],
            percentage_threshold=config["generate_suspectList"]["percentage_threshold"],
            min_depth=config["peak_calling"]["min_depth"]
        ) if config["peak_calling"]["filter_by_depth"] else  expand(
            "{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/round_{round}.SL.bed", 
            result_dir=config["results_dir"],
            num_reads=config["DownSample"]["Reads"],
            round=range(1, config["bootstrap_samples"]["rounds"] + 1),
            #round=range(1, 2 + 1),
            seacr_threshold=config["peak_calling"]["threshold"],
            percentage_threshold=config["generate_suspectList"]["percentage_threshold"]
            ),
       
       
        expand(
            "{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/all_rounds.SL.bed",
            result_dir=config["results_dir"],
            num_reads=config["DownSample"]["Reads"],
            seacr_threshold=config["peak_calling"]["threshold"],
            percentage_threshold=config["generate_suspectList"]["percentage_threshold"],
            min_depth=config["peak_calling"]["min_depth"]
        ) if config["peak_calling"]["filter_by_depth"] else  expand(
           "{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/all_rounds.SL.bed",
            result_dir=config["results_dir"],
            num_reads=config["DownSample"]["Reads"],
            seacr_threshold=config["peak_calling"]["threshold"],
            percentage_threshold=config["generate_suspectList"]["percentage_threshold"]
            ),


##############
# VALIDATION #
##############
# The following rules are for validating the suspect lists generated by the pipeline
#   call peaks for validation samples without removing SL regions
#   call peaks for validation samples after removing SL regions
#   overlap each validation sample peaks with SL generated by rule Initialize_SuspectPeak_Hunter

rule Validate_SL:
    input:
        # Call peaks for validation samples and overlap them to generate peaks in validation samples. It Generates peak files of validation samples without removing SL regions
        expand("{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/round_{round}.SL.bed", 
                result_dir=config["results_dir"],
                num_reads=config["DownSample"]["Reads"],
                round="Validation",
                seacr_threshold=config["peak_calling"]["threshold"],
                percentage_threshold=config["generate_suspectList"]["percentage_threshold"],
                min_depth=config["peak_calling"]["min_depth"]
                ) if config["peak_calling"]["filter_by_depth"] else  expand(
                "{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/round_{round}.SL.bed", 
                result_dir=config["results_dir"],
                num_reads=config["DownSample"]["Reads"],
                round="Validation",
                seacr_threshold=config["peak_calling"]["threshold"],
                percentage_threshold=config["generate_suspectList"]["percentage_threshold"]
                ),

        # Call Peaks for Validation Samples after removing SL
        expand("{result_dir}/05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}/noSL.{sample}.{seacr_mode}.bed",
                result_dir=config["results_dir"],
                round="Validation",
                num_reads=config["DownSample"]["Reads"],
                min_depth=config["peak_calling"]["min_depth"],
                seacr_threshold=config["peak_calling"]["threshold"],
                seacr_mode=config["peak_calling"]["mode"],
                sample=get_outfiles(test_samps=None,test_samps_st=None,type="all",array="validation_samples.tsv")
                ) if config["peak_calling"]["filter_by_depth"] else  expand(
                "{result_dir}/05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/SEACR_{seacr_threshold}/noSL.{sample}.{seacr_mode}.bed",
                result_dir=config["results_dir"],
                round="Validation",
                num_reads=config["DownSample"]["Reads"],
                min_depth=config["peak_calling"]["min_depth"],
                seacr_threshold=config["peak_calling"]["threshold"],
                seacr_mode=config["peak_calling"]["mode"],
                sample=get_outfiles(test_samps=None,test_samps_st=None,type="all",array="validation_samples.tsv")),

        # Overlap each Validation sample peaks with SL generated by rule Initialize_SuspectPeak_Hunter
        expand("{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/Overlap.Validation_Samples.SL.bed",
                result_dir=config["results_dir"],
                num_reads=config["DownSample"]["Reads"],
                min_depth=config["peak_calling"]["min_depth"],
                seacr_threshold=config["peak_calling"]["threshold"],
                percentage_threshold=config["generate_suspectList"]["percentage_threshold"]
                )if config["peak_calling"]["filter_by_depth"] else expand(
                "{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/Overlap.Validation_Samples.SL.bed",
                result_dir=config["results_dir"],
                num_reads=config["DownSample"]["Reads"],
                min_depth=config["peak_calling"]["min_depth"],
                seacr_threshold=config["peak_calling"]["threshold"],
                percentage_threshold=config["generate_suspectList"]["percentage_threshold"])





####################################################################
# (1) #  Prepare Genome for Mapping                                #
####################################################################


# genome index using bowtie-build to be used for bowtie mapping
rule genome_bowtie_index:
    input:
        reference="00.data/00.Genome/{genome}.fa"
    output:
        [
        "00.data/00.Genome/bowtie2_index/{genome}.1.bt2",
        "00.data/00.Genome/bowtie2_index/{genome}.2.bt2",
        "00.data/00.Genome/bowtie2_index/{genome}.3.bt2",
        "00.data/00.Genome/bowtie2_index/{genome}.4.bt2",
        "00.data/00.Genome/bowtie2_index/{genome}.rev.1.bt2",
        "00.data/00.Genome/bowtie2_index/{genome}.rev.2.bt2",
        ]
    conda:
        "envs/bowtie2.yaml"
    envmodules:
        "gcc",
        "bowtie2"
    params:
        path="00.data/00.Genome/bowtie2_index"
    log:
        "logs/02.00.{genome}_bowtie_index.log",
    threads: 32
    benchmark: "benchmarks/02.bowtie2_index.{genome}.tsv"
    shell:
        """
        bowtie2-build --threads {threads} {input.reference} {params.path}/{genome} 2> {log}
        """


# Genome fasta file index using samtools faidx
rule create_genome_fasta_index:
    input:
        reference="00.data/00.Genome/{genome}.fa"
    output:
        genome_index="00.data/00.Genome/{genome}.fa.fai"
    threads: 1
    envmodules:
        "gcc",
        "samtools"
    log: 
        "logs/00.Genome_Index.{genome}.log"
    benchmark: "benchmarks/00.Genome_Index.{genome}.tsv"
    shell:
        """
        samtools faidx {input.reference} -o {output.genome_index}
        """




####################################################################
# (2) #  FastQC and trimming                                       #
####################################################################

rule trim_galore_pe:
    input:
        [
            "00.data/00.rawReads/{sample}_1.fastq.gz", 
            "00.data/00.rawReads/{sample}_2.fastq.gz"
        ],
    output:
        fasta_fwd=temp("{result_dir}/01.trim/00.trimmedReads/{sample}_1.PE.fq.gz") 
        if config['keep_trimmed_reads_fastq'] != "True" 
        else "{result_dir}/01.trim/00.trimmedReads/{sample}_1.PE.fq.gz",
        
        fasta_rev=temp("{result_dir}/01.trim/00.trimmedReads/{sample}_2.PE.fq.gz") 
        if config['keep_trimmed_reads_fastq'] != "True" 
        else "{result_dir}/01.trim/00.trimmedReads/{sample}_2.PE.fq.gz",
        
        report_fwd="{result_dir}/01.trim/00.Reports/{sample}_1.PE.trimming_report.txt",
        report_rev="{result_dir}/01.trim/00.Reports/{sample}_2.PE.trimming_report.txt",
        
    threads: 48
    benchmark: "{result_dir}/benchmarks/00.trimming.{sample}.tsv"
    params:
        extra=lambda wildcards: check_if_adaptors_provided(
            st=all_samples_df, 
            sample=wildcards.sample, 
            trim_param=config["trim_galore_pe"]["additional_params"]
        )
    log:
        "{result_dir}/logs/01.trim_galore_{sample}.log",
    wrapper:
        "v3.9.0/bio/trim_galore/pe"

rule trim_galore_se:
    input:
        "00.data/00.rawReads/{sample}.fastq.gz",

    output:
        fasta=temp("{result_dir}/01.trim/00.trimmedReads/{sample}.SE.fq.gz") 
        if config['keep_trimmed_reads_fastq'] != "True" 
        else "{result_dir}/01.trim/00.trimmedReads/{sample}.SE.fq.gz",

        report="{result_dir}/01.trim/00.Reports/{sample}.SE.trimming_report.txt",

    threads: 48

    benchmark: "{result_dir}/benchmarks/00.trimming.{sample}.tsv"

    params:
        extra=config["trim_galore_se"]["additional_params"],

    log:
        "{result_dir}/logs/01.trim_galore_{sample}.log",

    wrapper:
        "v3.10.2/bio/trim_galore/se"


rule trimmedReads_fastqc_PE:
    input:
        trimmedread="{result_dir}/01.trim/00.trimmedReads/{sample}_{id}.PE.fq.gz"
    output:
        zip="{result_dir}/01.trim/01.fastqc/{sample}_{id}.PE_fastqc.zip",
        html="{result_dir}/01.trim/01.fastqc/{sample}_{id}.PE_fastqc.html"
    threads:
        24
    params:
        path="{result_dir}/01.trim/01.fastqc/"
    log:
        "{result_dir}/logs/01.trimmed_fastq_{sample}_{id}.log"
    benchmark: "{result_dir}/benchmarks/01.trimmed_fastq_{sample}_{id}.tsv"
    shell:
        """
        fastqc {input.trimmedread} --threads {threads} -o {params.path} 2> {log}
        """

rule trimmedReads_fastqc_SE:
    input:
        trimmedread="{result_dir}/01.trim/00.trimmedReads/{sample}.SE.fq.gz"
    output:
        zip="{result_dir}/01.trim/01.fastqc/{sample}.SE_fastqc.zip",
        html="{result_dir}/01.trim/01.fastqc/{sample}.SE_fastqc.html"
    threads:
        24
    params:
        path="{result_dir}/01.trim/01.fastqc/"
    log:
        "{result_dir}/logs/01.trimmed_fastq_{sample}.SE.log"
    benchmark: "{result_dir}/benchmarks/01.trimmed_fastq_{sample}.SE.tsv"
    shell:
        """
        fastqc {input.trimmedread} --threads {threads} -o {params.path} 2> {log}
        """




####################################################################
# (3) #  Mapping                                                   #
####################################################################
# Map to raw files if skip_trimming set to true in configfile else map to trimmed files

rule map_bowtie2_pe:
    input:
        read1="00.data/00.rawReads/{sample}_1.fastq.gz" if config["trim_galore_pe"]["skip_trimming"] == "True" else rules.trim_galore_pe.output.fasta_fwd,
        read2="00.data/00.rawReads/{sample}_2.fastq.gz" if config["trim_galore_pe"]["skip_trimming"] == "True" else rules.trim_galore_pe.output.fasta_rev,
        fasta=expand("00.data/00.Genome/bowtie2_index/{genome}.{id}.bt2",genome=config['genome'], id=["1","2","3","4","rev.1","rev.2"])
    output:
        bam=temp("{result_dir}/02.mapping/00.Raw/{sample}.sorted.bam") if config["keep_Raw_BAMs"] != "True" else "{result_dir}/02.mapping/00.Raw/{sample}.sorted.bam"
    conda:
        "envs/bowtie2.yaml"
    envmodules:
        "gcc",
        "bowtie2",
        "samtools"
    params:
        basename="{result_dir}/02.mapping/00.Raw/{sample}",
        index_path="00.data/00.Genome/bowtie2_index",
        genome="Danio_rerio.GRCz11.dna_sm.primary_assembly",
        additional_params=config["map_bowtie2_pe"]["additional_params"]
    threads: 48
    benchmark: "{result_dir}/benchmarks/02.bowtie2_mapping.{sample}.tsv"
    log:
        "{result_dir}/logs/02.01.mapping.{sample}.log",
    shell:
       """
       bowtie2 \
       {params.additional_params} \
       --threads={threads} \
       -x {params.index_path}/{params.genome} \
       -1 {input.read1} -2 {input.read2} |\
       samtools sort --threads {threads} -T {params.basename} -o {params.basename}.sorted.bam -
       """ 

rule map_bowtie2_se:
    input:
        read1="00.data/00.rawReads/{sample}.fastq.gz" if config["trim_galore_se"]["skip_trimming"] == "True" else "{result_dir}/01.trim/00.trimmedReads/{sample}.SE.fq.gz",
        fasta=expand("00.data/00.Genome/bowtie2_index/{genome}.{id}.bt2",genome=config['genome'], id=["1","2","3","4","rev.1","rev.2"])
    output:
        bam=temp("{result_dir}/02.mapping/00.Raw/{sample}.sorted.bam") if config["keep_Raw_BAMs"] != "True" else "{result_dir}/02.mapping/00.Raw/{sample}.sorted.bam"
    conda:
        "envs/bowtie2.yaml"
    envmodules:
        "gcc",
        "bowtie2",
        "samtools"
    params:
        basename="{result_dir}/02.mapping/00.Raw/{sample}",
        index_path="00.data/00.Genome/bowtie2_index",
        genome="Danio_rerio.GRCz11.dna_sm.primary_assembly",
        additional_params=config["map_bowtie2_se"]["additional_params"]
    threads: 48
    benchmark: "{result_dir}/benchmarks/02.bowtie2_mapping.{sample}.tsv"
    log:
        "{result_dir}/logs/02.01.mapping.{sample}.log"
    shell:
        """
        bowtie2 \
        {params.additional_params} \
        --threads={threads} \
        -U {input.read1} \
        -x {params.index_path}/{params.genome} |\
        samtools sort \
        --threads {threads} \
        -T {params.basename} \
        -o {params.basename}.sorted.bam -
        """

####################################################################
# (4) #  Remove mitochondrial chromosome                           #
####################################################################

rule remove_mitochondrial_chromosome:
    input:
        bam="{result_dir}/02.mapping/00.Raw/{sample}.sorted.bam"
    output:
        intermediate=temp(directory("{result_dir}/02.mapping/00.Raw/{sample}")),
        bam_wo_MT=temp("{result_dir}/02.mapping/00.Raw/{sample}.no_MT.sorted.bam") if config["keep_BAMs_wo_MT"] != "True" else "{result_dir}/02.mapping/00.Raw/{sample}.no_MT.sorted.bam"
    threads: 8
    envmodules:
        "gcc",
        "samtools"
    log:
        "{result_dir}/logs/02.03.remove_MT.{sample}.log"
    benchmark: "{result_dir}/benchmarks/02.03.remove_MT.{sample}.tsv"
    message:
        "**********Running Rule: remove_mitochondrial_chromosome************"
    shell:
        "scripts/Remove_Mitochondrial_Chromosome.sh {wildcards.sample} {input.bam} {output.intermediate} {output.bam_wo_MT} {threads}"







## TO BE ADDED ##  #### COUNT READS LEFT AFTER REMOVING MITOCHONDRIAL READS <-- done in the "remove_mitochondrial_chromosome" rule

####################################################################
# (5) #  BootStrapping                                             #
####################################################################


#*******#
# (5.1) #  Downsample samples for each bootstrap round
#*******#


rule downsample_bam:
    input:
        raw="{result_dir}/02.mapping/00.Raw/{sample}.no_MT.sorted.bam",
        #bam="05.bootstrapping/round_{round}/00.mapping/{sample}.no_MT.sorted.bam"
        bam=rules.remove_mitochondrial_chromosome.output.bam_wo_MT,
        #sample_array="04.bootstrapping_arrays/array_{round}.tsv",
    output:
        #downsampled_bam="02.mapping/01.DownSampling.{num_reads}/{round}/{sample}.no_MT.sorted.bam"
        downsampled_bam=temp("{result_dir}/05.bootstrapping/round_{round}/01.DownSampling.{num_reads}/{sample}.no_MT.sorted.bam") if config["keep_downsampled_bam"] != "True" else "{result_dir}/05.bootstrapping/round_{round}/01.DownSampling.{num_reads}/{sample}.no_MT.sorted.bam"
    conda:
        "envs/bowtie2.yaml"
    envmodules:
        "gcc",
        "samtools"
    params:
        reads=config["DownSample"]["Reads"],
        reads_PE=config["DownSample"]["Reads_PE"],
        reads_SE=config["DownSample"]["Reads_SE"],
        sample_type=lambda wildcards: get_type(wildcards.sample, all_samples_df)
    threads: config["DownSample"]["threads"]
    benchmark:
        #"benchmarks/02.02.downsample_bam.round_{round}.{sample}.tsv"
        "{result_dir}/benchmarks/02.02.downsample_bam.{num_reads}.round_{round}.{sample}.tsv"
    log:
        #"logs/02.02.downsample_bam.{round}.{sample}.log"
        "{result_dir}/logs/02.02.downsample_bam.{num_reads}.round_{round}.{sample}.log"
    shell:
        "scripts/downsample_bams.sh {input.bam} {params.reads_PE} {params.reads_SE} {params.sample_type} {output.downsampled_bam} {threads}"

#*******#
# (5.2) #  BAM2BED 
#*******#


rule bamtobed:
    input: 
        #bam="02.mapping/01.DownSampling.{num_reads}/{sample}.no_MT.sorted.bam"
        bam=rules.downsample_bam.output.downsampled_bam
    output:
        #intermediate=temp(directory("03.PeakCalling/00.rawData/DownSampling_{num_reads}/intermediate/{sample}")),
        #bamtobed=temp("03.PeakCalling/00.rawData/DownSampling_{num_reads}/00.bamtobed/{sample}.bamtobed")
        intermediate=temp(directory("{result_dir}/05.bootstrapping/round_{round}/02.bamtobed.DS_{num_reads}/intermediate/{sample}")),
        bamtobed=temp("{result_dir}/05.bootstrapping/round_{round}/02.bamtobed.DS_{num_reads}/{sample}.bamtobed")
    threads: 32
    params:
        sample_type=lambda wildcards: get_type(wildcards.sample, all_samples_df)
    envmodules:
        "gcc",
        "bedtools2",
        "samtools"
    log: 
       "{result_dir}/logs/03.00.round_{round}.Downsampling_{num_reads}.bamtobed.{sample}.log"
    benchmark: "{result_dir}/benchmarks/03.00.round_{round}.Downsampling_{num_reads}.bamtobed.{sample}.tsv"
    shell:
        "scripts/04.01.GenomeCov_PE_SE.sh {wildcards.sample} {input.bam} {output.intermediate} {output.bamtobed} {threads} {params.sample_type} > {log} 2>&1"


#*******#
# (5.3) #  GenomeCov
#*******#

rule genome_cov:
    input: 
        bamtobed=rules.bamtobed.output.bamtobed,
    output:
        genomecov_out=temp(
            "{result_dir}/05.bootstrapping/round_{round}/03.genomecov.DS_{num_reads}/{sample}.bedgraph"
        ) if config["keep_genomecov_bedgraph"] != "True" else 
            "{result_dir}/05.bootstrapping/round_{round}/03.genomecov.DS_{num_reads}/{sample}.bedgraph",
    threads: 8
    log: "{result_dir}/logs/03.01.round_{round}.Downsampling_{num_reads}.genome_cov.{sample}.log"
    benchmark: "{result_dir}/benchmarks/03.01.round_{round}.Downsampling_{num_reads}.genome_cov.{sample}.tsv"
    params:
        genome_name=genome,
        genome_path="00.data/00.Genome",
    envmodules:
        "gcc",
        "bedtools2"
    shell:
        """
        bedtools genomecov -bg -i {input.bamtobed} -g {params.genome_path}/{params.genome_name}.fa.fai > {output.genomecov_out}
        """



#*****************#
# (5.4* OPTIONAL) #  Filter GenomeCov (Optional)
#*****************#



rule filter_genome_cov:
    input: 
        genomecov=rules.genome_cov.output.genomecov_out,
        
    output:
        genomecov_out_filtered=temp(
            "{result_dir}/05.bootstrapping/round_{round}/03.genomecov.DS_{num_reads}/{sample}.filtered_mindepth_{min_depth}.bedgraph"
        ) if config["keep_genomecov_bedgraph"] != "True" else 
            "{result_dir}/05.bootstrapping/round_{round}/03.genomecov.DS_{num_reads}/{sample}.filtered_mindepth_{min_depth}.bedgraph",
        
    threads: 8
    log: "{result_dir}/logs/03.01.round_{round}.Downsampling_{num_reads}.genome_cov.filtered_mindepth_{min_depth}.{sample}.log"
    benchmark: "{result_dir}/benchmarks/03.01.round_{round}.Downsampling_{num_reads}.genome_cov.filtered_mindepth_{min_depth}.{sample}.tsv"
    params:
          genome_name=genome,
          genome_path="00.data/00.Genome",
          min_depth=config["peak_calling"]["min_depth"],
    envmodules:
        "gcc",
        "bedtools2"
    shell:
        """
        awk -v min_depth={params.min_depth} '{{if($4 > min_depth) print}}' {input.genomecov} > {output.genomecov_out_filtered}
        """


#*******#
# (5.5) #  PeakCalling
#*******#


seacr_mode=config["peak_calling"]["mode"]
rule peak_calling:
    input:
        bebdgraph_file=rules.filter_genome_cov.output.genomecov_out_filtered 
        if config["peak_calling"]["filter_by_depth"] else 
        rules.genome_cov.output.genomecov_out
    output:
        seacr_out="{result_dir}/05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}/{sample}.{seacr_mode}.bed" 
        if config["peak_calling"]["filter_by_depth"] else 
        "{result_dir}/05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/SEACR_{seacr_threshold}/{sample}.{seacr_mode}.bed",
    params:
        threshold=config["peak_calling"]["threshold"],
        normalization=config["peak_calling"]["normalization"],
        mode=config["peak_calling"]["mode"],
        output_prefix="{result_dir}/05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}/{sample}" 
        if config["peak_calling"]["filter_by_depth"] else 
        "{result_dir}/05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/SEACR_{seacr_threshold}/{sample}",
    log:  
       "{result_dir}/logs/03.02.round_{round}.Downsampling_{num_reads}.peak_calling.{sample}.filtered_mindepth_{min_depth}.{seacr_mode}.SEACR_{seacr_threshold}.log" 
       if config["peak_calling"]["filter_by_depth"] else 
       "{result_dir}/logs/03.02.round_{round}.Downsampling_{num_reads}.peak_calling.{sample}.{seacr_mode}.SEACR_{seacr_threshold}.log"

    shell:
       """
       echo "---- {wildcards.sample}:00:Peak Calling using threshold: {params.threshold} normalization:{params.normalization} mode:{params.mode} ------"
       SEACR_1.3.sh {input.bebdgraph_file} {params.threshold} {params.normalization} {params.mode} {params.output_prefix}  > {log} 2>&1
       """

#**************************#
# (5.6) #  Overlap peaks   #
#**************************#


rule overlap_peaks:
    input:
        #sample_array="04.bootstrapping_arrays/array_{round}.tsv",
        sample_array=lambda wildcards: "04.bootstrapping_arrays/array_{}.tsv".format(wildcards.round) 
                 if wildcards.round != "Validation" 
                 else "validation_samples.tsv",
        #expanded_files=rules.bootstrap_samples.output.expanded_files,
        peaks=lambda wildcards: expand("{result_dir}/05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/SEACR_{seacr_threshold}/{sample}.{seacr_mode}.bed", 
                                        results_dir=wildcards.results_dir,
                                        round=wildcards.round, 
                                        seacr_threshold=wildcards.seacr_threshold,
                                        num_reads=wildcards.num_reads, 
                                        sample=get_outfiles(test_samps=None, 
                                                            test_samps_st=None, 
                                                            type="all",
                                                            array="04.bootstrapping_arrays/array_{}.tsv".format(wildcards.round) 
                                                            if wildcards.round != "Validation" 
                                                            else "validation_samples.tsv"),
                                        seacr_mode=config["peak_calling"]["mode"]),
    output:
        merged_regions_binary_data = "{result_dir}/05.bootstrapping/05.Overlap_Peaks.DS_{num_reads}/SEACR_{seacr_threshold}.mergepeaks/round_{round}.rawbam.mergedPeaks.binary.tab",
        merged_regions_collapsed = "{result_dir}/05.bootstrapping/05.Overlap_Peaks.DS_{num_reads}/SEACR_{seacr_threshold}.mergepeaks/round_{round}.rawbam.mergedPeaks.tab"
    threads: 1
    params:
        #peak_dir="03.PeakCalling/00.rawData/DownSampling_{num_reads}/02.peaks",
        peak_dir="{result_dir}/05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/SEACR_{seacr_threshold}",
        bedfile_suffix=expand("{mode}.bed",mode=config["peak_calling"]["mode"]),
        #sample_array="04.bootstrapping_arrays/array_{round}.tsv"
    envmodules:
        "gcc",
        "bedtools2"
    log:
        "{result_dir}/logs/04.00.round_{round}.DownSampling_{num_reads}.SEACR_{seacr_threshold}.round_{round}.overlap_peaks.log"
    benchmark: "{result_dir}/benchmarks/04.00.round_{round}.DownSampling_{num_reads}.SEACR_{seacr_threshold}.round_{round}.overlap_peaks.tsv"
    shell:
        "scripts/MergeSelectedPeakFiles.sh {input.sample_array} {params.peak_dir} {params.bedfile_suffix} {output.merged_regions_binary_data} {output.merged_regions_collapsed}"




use rule overlap_peaks as overlap_peaks_from_filtered_genome_cov with:
    input:
        #sample_array="04.bootstrapping_arrays/array_{round}.tsv" if "{round}"!="Validation" else "validation_samples.tsv",
        sample_array=lambda wildcards: "04.bootstrapping_arrays/array_{}.tsv".format(wildcards.round) 
                 if wildcards.round != "Validation" 
                 else "validation_samples.tsv",
        #expanded_files=rules.bootstrap_samples.output.expanded_files,
        peaks=lambda wildcards: expand("{result_dir}/05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}/{sample}.{seacr_mode}.bed", 
                                        result_dir=wildcards.result_dir,
                                        round=wildcards.round, 
                                        seacr_threshold=wildcards.seacr_threshold,
                                        num_reads=wildcards.num_reads, 
                                        #sample=get_outfiles(test_samps=None, test_samps_st=None, type="all",array="04.bootstrapping_arrays/array_{}.tsv".format(wildcards.round)),
                                        sample=get_outfiles(test_samps=None, 
                                                            test_samps_st=None, 
                                                            type="all",
                                                            array="04.bootstrapping_arrays/array_{}.tsv".format(wildcards.round) 
                                                            if wildcards.round != "Validation" 
                                                            else "validation_samples.tsv"),
                                        seacr_mode=config["peak_calling"]["mode"],
                                        min_depth=wildcards.min_depth),
    output:
        merged_regions_binary_data = "{result_dir}/05.bootstrapping/05.Overlap_Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.mergepeaks/round_{round}.rawbam.mergedPeaks.binary.tab",
        merged_regions_collapsed = "{result_dir}/05.bootstrapping/05.Overlap_Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.mergepeaks/round_{round}.rawbam.mergedPeaks.tab"
    params:
        peak_dir="{result_dir}/05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}",
        bedfile_suffix=expand("{mode}.bed",mode=config["peak_calling"]["mode"]),
        #sample_array="04.bootstrapping_arrays/array_{round}.tsv"
    log:
        "{result_dir}/logs/04.00.round_{round}.DownSampling_{num_reads}.filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.round_{round}.overlap_peaks.log"
    benchmark: "{result_dir}/benchmarks/04.00.round_{round}.DownSampling_{num_reads}.filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.round_{round}.overlap_peaks.tsv"




#************************************#
# (5.7) #  Generating Suspect Lists  #
#************************************#



percentage_threshold=config["generate_suspectList"]["percentage_threshold"]
rule generate_suspectList:
    input:
       raw_peaks= rules.overlap_peaks_from_filtered_genome_cov.output.merged_regions_binary_data if config["peak_calling"]["filter_by_depth"] else rules.overlap_peaks.output.merged_regions_binary_data,
       chr_lengths="00.data/00.Genome/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.fai",
    output:
        SuspectList="{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/round_{round}.SL.bed" 
        if config["peak_calling"]["filter_by_depth"] else 
        "{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/round_{round}.SL.bed",
    threads: 1
    params:
        num_of_samples=config["generate_suspectList"]["num_of_samples"],
        minimum_length=config["generate_suspectList"]["minimum_length"],
        percentage_threshold=config["generate_suspectList"]["percentage_threshold"],
        metadata=SL_generation_samples,
        outpath="{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}" 
        if config["peak_calling"]["filter_by_depth"] else 
        "{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}"
    envmodules:
        "gcc",
        "bedtools2",
        "r"
    log:
       "{result_dir}/logs/06.Generating_SuspectLists.DS_{num_reads}.filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}.round_{round}.log" 
       if config["peak_calling"]["filter_by_depth"] else 
       "{result_dir}/logs/06.Generating_SuspectLists.DS_{num_reads}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}.round_{round}.log"
    shell:
       """
       module load gcc bedtools2
       mkdir -p {params.outpath}
       Rscript scripts/Generate_SuspectLists.bootstrap.R \
       --raw_peaks_file={input.raw_peaks} \
       --chr_lengths_file={input.chr_lengths} \
       --metadata={params.metadata} \
       --num_sample={params.num_of_samples} \
       --percentage_threshold={params.percentage_threshold} \
       --min_shared_region_len={params.minimum_length} \
       --output_BL_bases={params.outpath}/round_{wildcards.round}.SL \
       --ignore_filtered_plots="TRUE" > {log} 2>&1

       sort -k1,1V -k2,2n {output.SuspectList} > {params.outpath}/round_{wildcards.round}.SL.sorting.bed
       mv {params.outpath}/round_{wildcards.round}.SL.sorting.bed {output.SuspectList}

       """

# BOOTSTRAPPING ENDS HERE

######################################
# (6) #  Overlap Suspect Lists       #
######################################

rule Regions_Present_in_all_SL:
    input:
       SL= expand("{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/round_{round}.SL.bed", 
        result_dir=config["results_dir"],
        num_reads=config["DownSample"]["Reads"],
        round=range(1, config["bootstrap_samples"]["rounds"] + 1),
        #round=range(1, 2 + 1),
        seacr_threshold=config["peak_calling"]["threshold"],
        percentage_threshold=config["generate_suspectList"]["percentage_threshold"],
        min_depth=config["peak_calling"]["min_depth"])
        if config["peak_calling"]["filter_by_depth"] else
        expand("{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/round_{round}.SL.bed", 
        result_dir=config["results_dir"],
        num_reads=config["DownSample"]["Reads"],
        round=range(1, config["bootstrap_samples"]["rounds"] + 1),
        #round=range(1, 2 + 1),
        seacr_threshold=config["peak_calling"]["threshold"],
        percentage_threshold=config["generate_suspectList"]["percentage_threshold"]),
       #filt_peaks=rules.overlap_peaks_filtered.output.merged_regions_binary_data,
       chr_lengths="00.data/00.Genome/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.fai",
    output:
        Overlap_SLs="{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/all_rounds.SL.bed" 
        if config["peak_calling"]["filter_by_depth"] else 
        "{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/all_rounds.SL.bed",
        
        Overlap_SLs_filtered="{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/all_rounds.SL.filtered.bed" 
        if config["peak_calling"]["filter_by_depth"] else 
        "{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/all_rounds.SL.filtered.bed"
    params:
        input_dir="{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}" 
        if config["peak_calling"]["filter_by_depth"] else 
        "{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}",
        suspectlist_min_bootstraps=config["generate_suspectList"]["min_num_of_bootstrap_rounds"]
    threads: 5
    envmodules:
        "gcc",
        "bedtools2",
        "r"
    log:
       "{result_dir}/logs/05.00.DownSampled_{num_reads}.filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}.log"
       if config["peak_calling"]["filter_by_depth"] else
       "{result_dir}/logs/05.00.DownSampled_{num_reads}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}.log"
    shell:
       """
        scripts/overlap_SL.sh {params.input_dir} {params.suspectlist_min_bootstraps} {output.Overlap_SLs} {output.Overlap_SLs_filtered}  
       """


#########################################
### (7) ValidateSL                    ###
#########################################



#*******#
# (7.1) #  Remove SuspectList regions from BAMTOBED files
#*******#

rule remove_Suspect_Regions_from_BAMTOBED:
   input: 
      bamtobed=rules.bamtobed.output.bamtobed,
      #suspect_list=rules.Regions_Present_in_all_SL.output.Overlap_SLs_filtered
      suspect_list=expand(
            "{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/all_rounds.SL.filtered.bed",
            result_dir=config["results_dir"],
            num_reads=config["DownSample"]["Reads"],
            seacr_threshold=config["peak_calling"]["threshold"],
            percentage_threshold=config["generate_suspectList"]["percentage_threshold"],
            min_depth=config["peak_calling"]["min_depth"]
        ) if config["peak_calling"]["filter_by_depth"] else  expand(
           "{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/all_rounds.SL.filtered.bed",
            result_dir=config["results_dir"],
            num_reads=config["DownSample"]["Reads"],
            seacr_threshold=config["peak_calling"]["threshold"],
            percentage_threshold=config["generate_suspectList"]["percentage_threshold"]
            ),
   output: 
      bamtobed_filtered=temp("{result_dir}/05.bootstrapping/round_{round}/02.bamtobed.DS_{num_reads}/noSL.{sample}.bamtobed")
   threads: 12

   log:   "{result_dir}/logs/03.00.round_{round}.Downsampling_{num_reads}.bamtobed.{sample}.noSL.log"
   shell: 
      """
      echo "---- Removing suspect regions ------"
      bedtools intersect -a {input.bamtobed} -b {input.suspect_list} -v > {output.bamtobed_filtered} 
      """


#*******#
# (7.2) #  Calculate genome covarage
#*******#

rule genomecov_SL_Removed:
    input: 
        #bamtobed="03.PeakCalling/00.rawData/DownSampling_{num_reads}/00.bamtobed/{sample}.bamtobed"
        bamtobed=rules.remove_Suspect_Regions_from_BAMTOBED.output.bamtobed_filtered
    output:
        genomecov_out=temp(
            "{result_dir}/05.bootstrapping/round_{round}/03.genomecov.DS_{num_reads}/noSL.{sample}.bedgraph"
        ) if config["keep_genomecov_bedgraph"] != "True" else 
            "{result_dir}/05.bootstrapping/round_{round}/03.genomecov.DS_{num_reads}/noSL.{sample}.bedgraph",
    threads: 8
    log: "{result_dir}/logs/03.01.round_{round}.Downsampling_{num_reads}.genome_cov.{sample}.noSL.log"
    benchmark: "{result_dir}/benchmarks/03.01.round_{round}.Downsampling_{num_reads}.genome_cov.{sample}.noSL.tsv"
    params:
          genome_name=genome,
          genome_path="00.data/00.Genome"
    envmodules:
        "gcc",
        "bedtools2"
    shell:
        """
        bedtools genomecov -bg -i {input.bamtobed} -g {params.genome_path}/{params.genome_name}.fa.fai > {output.genomecov_out}
        """


#*******#
# (7.3) #  Filter Genomecov
#*******#
rule filter_genome_cov_SL_Removed:
    input: 
        genomecov=rules.genomecov_SL_Removed.output.genomecov_out
    output:
        genomecov_out_filtered=temp(
            "{result_dir}/05.bootstrapping/round_{round}/03.genomecov.DS_{num_reads}/noSL.{sample}.filtered_mindepth_{min_depth}.bedgraph"
        ) if config["keep_genomecov_bedgraph"] != "True" else 
            "{result_dir}/05.bootstrapping/round_{round}/03.genomecov.DS_{num_reads}/noSL.{sample}.filtered_mindepth_{min_depth}.bedgraph"

    threads: 8
    log: "{result_dir}/logs/03.01.round_{round}.Downsampling_{num_reads}.genome_cov.filtered_mindepth_{min_depth}.{sample}.noSL.log"
    benchmark: "{result_dir}/benchmarks/03.01.round_{round}.Downsampling_{num_reads}.genome_cov.filtered_mindepth_{min_depth}.{sample}.noSL.tsv"
    params:
          genome_name=genome,
          genome_path="00.data/00.Genome",
          min_depth=config["peak_calling"]["min_depth"]
    envmodules:
        "gcc",
        "bedtools2"
    shell:
        """
        awk -v min_depth={params.min_depth} '{{if($4 > min_depth) print}}' {input.genomecov} > {output.genomecov_out_filtered}
        """

#*******#
# (7.4) #  Call Peaks
#*******#
rule peak_calling_SL_Removed:
    input:
        bebdgraph_file=rules.filter_genome_cov_SL_Removed.output.genomecov_out_filtered 
        if config["peak_calling"]["filter_by_depth"] else 
        rules.genomecov_SL_Removed.output.genomecov_out
    output:
        seacr_out="{result_dir}/05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}/noSL.{sample}.{seacr_mode}.bed" 
        if config["peak_calling"]["filter_by_depth"] else 
        "{result_dir}/05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/SEACR_{seacr_threshold}/noSL.{sample}.{seacr_mode}.bed",
    params:
        threshold=config["peak_calling"]["threshold"],
        normalization=config["peak_calling"]["normalization"],
        mode=config["peak_calling"]["mode"],
        output_prefix="{result_dir}/05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}/noSL.{sample}" 
        if config["peak_calling"]["filter_by_depth"] else 
        "{result_dir}/05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/SEACR_{seacr_threshold}/noSL.{sample}",
    log:  
       "{result_dir}/logs/03.02.round_{round}.Downsampling_{num_reads}.peak_calling.{sample}.filtered_mindepth_{min_depth}.{seacr_mode}.SEACR_{seacr_threshold}.noSL.log" 
       if config["peak_calling"]["filter_by_depth"] else 
       "{result_dir}/logs/03.02.round_{round}.Downsampling_{num_reads}.peak_calling.{sample}.{seacr_mode}.SEACR_{seacr_threshold}.{peak_filter}.log"
    shell:
       """
       echo "---- {wildcards.sample}:00:Peak Calling using threshold: {params.threshold} normalization:{params.normalization} mode:{params.mode} ------"
       SEACR_1.3.sh {input.bebdgraph_file} {params.threshold} {params.normalization} {params.mode} {params.output_prefix}  > {log} 2>&1
       """

#*******#
# (7.5) #  Overlap each validation smaple with SL
#*******#

rule Overlap_SLs_and_Validation_Samples:
    input:
        validation_samples="validation_samples.tsv",
        peaks=lambda wildcards: expand(
                "{result_dir}/05.bootstrapping/round_Validation/04.Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}/{sample}.{seacr_mode}.bed",
                result_dir=config["results_dir"],
                num_reads=config["DownSample"]["Reads"],
                seacr_threshold=config["peak_calling"]["threshold"],
                min_depth=config["peak_calling"]["min_depth"],
                sample=get_outfiles(test_samps=None, 
                                    test_samps_st=None, 
                                    type="all",
                                    array="validation_samples.tsv"),
                seacr_mode=config["peak_calling"]["mode"]
            ),
        bams=lambda wildcards: expand(
                "{result_dir}/05.bootstrapping/round_Validation/01.DownSampling.{num_reads}/{sample}.no_MT.sorted.bam",
                result_dir=config["results_dir"],
                num_reads=config["DownSample"]["Reads"],
                sample=get_outfiles(test_samps=None, 
                                    test_samps_st=None, 
                                    type="all",
                                    array="validation_samples.tsv")
            ),
        SL=rules.remove_Suspect_Regions_from_BAMTOBED.input.suspect_list
    output:
        overlap_VL_SL="{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/Overlap.Validation_Samples.SL.bed" 
        if config["peak_calling"]["filter_by_depth"] else 
        "{result_dir}/05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/Overlap.Validation_Samples.SL.bed",
    params:
        results_dir=config["results_dir"],
        num_reads=config["DownSample"]["Reads"],
        seacr_threshold=config["peak_calling"]["threshold"],
        min_depth=config["peak_calling"]["min_depth"],
        seacr_mode=config["peak_calling"]["mode"],
        #BED_FILES_SUFFIX="stringent.bed"
    shell:
        """

        SAMPLE_FILE={input.validation_samples}
        BED_FILE_DIR="{params.results_dir}/05.bootstrapping/round_Validation/04.Peaks.DS_{params.num_reads}/filtered_mindepth_{params.min_depth}.SEACR_{params.seacr_threshold}"
        SAMPLES=$(awk 'NR > 1 {{print $4 ".{params.seacr_mode}.bed"}}' {input.validation_samples} | sort | uniq)
        BED_FILES=$(for SAMPLE in $SAMPLES; do echo "$BED_FILE_DIR/$SAMPLE"; done)
        awk '{{filename=FILENAME; sub(/^.*\//, "", filename); print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" filename}}' $BED_FILES | sort -k1,1V -k2,2n | intersectBed -b stdin -a {input.SL} -wao > {output.overlap_VL_SL}

        """


