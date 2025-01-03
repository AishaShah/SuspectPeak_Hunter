# path to configuration file for input parameters
configfile: 'config.yaml'
workdir: config["workdir"]
# load required python packages
import pandas as pd

#***********************************************#
#          Helper Functions for Metadata        #
#***********************************************#

from helpers import *

#*************#
# INPUT FILES #
#*************#
genome = config["genome"] # Load genome


#*****#
# (1) #  Input files for generating suspect list
#*****#

# Load samplesheet for samples needed for suspectlist generation
SL_array=config["SL_array"] 
# read samplesheet and set "sample" name column as row index
st = pd.read_table(SL_array).set_index('sample', drop=False) 
samps = get_samples(st)   # get sample names
seacr_mode=config["peak_calling"]["mode"]
seacr_threshold=config["peak_calling"]["threshold"]



#*****#
# (2) #  Input files for applying suspect list
#*****#


SE=get_outfiles(samps, st, type="SE")
PE=get_outfiles(samps, st, type="PE")

##***************************************************##
## Generate SuspectList                              ##
##***************************************************##


# rule Initialize_SuspectPeak_Hunter runs all the other rules required for generating suspectlist
rule Initialize_SuspectPeak_Hunter:
    input:
        #expand("00.data/00.Genome/{genome}.fa.fai" , genome=config["genome"]),
        #expand("02.mapping/00.Raw/{sample}.no_MT.sorted.bam",sample=PE) ,
        #expand("02.mapping/00.Raw/{sample}.no_MT.sorted.bam",sample=SE),
        # overlap peaks per round
        expand(
            "05.bootstrapping/05.Overlap_Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.mergepeaks/round_{round}.rawbam.mergedPeaks.binary.tab", 
            num_reads=config["DownSample"]["Reads"],
            #round=range(1, config["bootstrap_samples"]["rounds"] + 1),
            round=range(1, 2 + 1),
            seacr_threshold=config["peak_calling"]["threshold"],
            min_depth=config["peak_calling"]["min_depth"]
         ) if config["peak_calling"]["filter_by_depth"] else  expand(
            "05.bootstrapping/05.Overlap_Peaks.DS_{num_reads}/SEACR_{seacr_threshold}.mergepeaks/round_{round}.rawbam.mergedPeaks.binary.tab", 
            num_reads=config["DownSample"]["Reads"],
            #round=range(1, config["bootstrap_samples"]["rounds"] + 1),
            round=range(1, 2 + 1),
            seacr_threshold=config["peak_calling"]["threshold"]
            ),
        
        # generate suspect list per round
        expand("05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/round_{round}.SL.bed", 
            num_reads=config["DownSample"]["Reads"],
            #round=range(1, config["bootstrap_samples"]["rounds"] + 1),
            round=range(1, 2 + 1),
            seacr_threshold=config["peak_calling"]["threshold"],
            percentage_threshold=config["generate_suspectList"]["percentage_threshold"],
            min_depth=config["peak_calling"]["min_depth"]
        ) if config["peak_calling"]["filter_by_depth"] else  expand(
            "05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/round_{round}.SL.bed", 
            num_reads=config["DownSample"]["Reads"],
            #round=range(1, config["bootstrap_samples"]["rounds"] + 1),
            round=range(1, 2 + 1),
            seacr_threshold=config["peak_calling"]["threshold"],
            percentage_threshold=config["generate_suspectList"]["percentage_threshold"]
            ),
       
       
        expand(
            "05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/all_rounds.SL.bed",
            num_reads=config["DownSample"]["Reads"],
            seacr_threshold=config["peak_calling"]["threshold"],
            percentage_threshold=config["generate_suspectList"]["percentage_threshold"],
            min_depth=config["peak_calling"]["min_depth"]
        ) if config["peak_calling"]["filter_by_depth"] else  expand(
           "05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/all_rounds.SL.bed",
            num_reads=config["DownSample"]["Reads"],
            seacr_threshold=config["peak_calling"]["threshold"],
            percentage_threshold=config["generate_suspectList"]["percentage_threshold"]
            ),




####################################################################
# (1) #  Setup Input Directory <-- create softlinks to inuput data #
####################################################################

# softlinks to files stored in data strage i.e nas cannot be accessed in slurm jobs 
# it is better to prepare these folders outside snakemake
# copy files or create soft links (if data stored on accessible device)
# rule setup_input_sample_directory:
#     input:
#         array=SL_array
#     output:
#         expand("00.data/00.rawReads/{sample}_{id}.fastq.gz",sample=PE, id=["1","2"]) ,
#         expand("00.data/00.rawReads/{sample}.fastq.gz",sample=SE) ,
#         directory("00.data/00.rawReads")
#     params:
#         fq_dir="00.data/00.rawReads",
#         fq_trimmed_dir="01.trim/00.trimmedReads",
#         bam_dir="02.mapping/00.Raw"
#     log:
#         "logs/00.setup_input_sample_directory.log"
#     shell:
#         """
#         scripts/Prepare_Input_Sample_directory.sh {input.array} {params.fq_dir} {params.fq_trimmed_dir} {params.bam_dir}
#         """



####################################################################
# (2) #  Prepare Genome for Mapping                                #
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
# (3) #  FastQC and trimming                                       #
####################################################################


rule trim_galore_pe:
    input:
        ["00.data/00.rawReads/{sample}_1.fastq.gz", "00.data/00.rawReads/{sample}_2.fastq.gz"],
    output:
        fasta_fwd="01.trim/00.trimmedReads/{sample}_1.PE.fq.gz",
        report_fwd="01.trim/00.Reports/{sample}_1.PE.trimming_report.txt",
        fasta_rev="01.trim/00.trimmedReads/{sample}_2.PE.fq.gz",
        report_rev="01.trim/00.Reports/{sample}_2.PE.trimming_report.txt",
    threads: 48
    benchmark: "benchmarks/00.trimming.{sample}.tsv"
    params:
        #extra=config["trim_galore_pe"]["additional_params"]
        extra=lambda wildcards: check_if_adaptors_provided(st=st, sample=wildcards.sample, trim_param=config["trim_galore_pe"]["additional_params"])
    log:
        "logs/01.trim_galore_{sample}.log",
    wrapper:
        "v3.9.0/bio/trim_galore/pe"

rule trim_galore_se:
    input:
        "00.data/00.rawReads/{sample}.fastq.gz",
    output:
        fasta="01.trim/00.trimmedReads/{sample}.SE.fq.gz",
        report="01.trim/00.Reports/{sample}.SE.trimming_report.txt",
    threads: 48
    benchmark: "benchmarks/00.trimming.{sample}.tsv"
    params:
        extra=config["trim_galore_se"]["additional_params"],
    log:
        "logs/01.trim_galore_{sample}.log",
    wrapper:
        "v3.10.2/bio/trim_galore/se"



## TO BE ADDED ## COUNT READS LEFT AFTER TRIMMING



rule trimmedReads_fastqc_PE:
    input:
        trimmedread="01.trim/00.trimmedReads/{sample}_{id}.PE.fq.gz"
    output:
        zip="01.trim/01.fastqc/{sample}_{id}.PE_fastqc.zip",
        html="01.trim/01.fastqc/{sample}_{id}.PE_fastqc.html"
    threads:
        24
    params:
        path="01.trim/01.fastqc/"
    log:
        "logs/01.trimmed_fastq_{sample}_{id}.log"
    benchmark: "benchmarks/01.trimmed_fastq_{sample}_{id}.tsv"
    shell:
        """
        fastqc {input.trimmedread} --threads {threads} -o {params.path} 2> {log}
        """

rule trimmedReads_fastqc_SE:
    input:
        trimmedread="01.trim/00.trimmedReads/{sample}.SE.fq.gz"
    output:
        zip="01.trim/01.fastqc/{sample}.SE_fastqc.zip",
        html="01.trim/01.fastqc/{sample}.SE_fastqc.html"
    threads:
        24
    params:
        path="01.trim/01.fastqc/"
    log:
        "logs/01.trimmed_fastq_{sample}.SE.log"
    benchmark: "benchmarks/01.trimmed_fastq_{sample}.SE.tsv"
    shell:
        """
        fastqc {input.trimmedread} --threads {threads} -o {params.path} 2> {log}
        """




####################################################################
# (4) #  Mapping                                                   #
####################################################################
# Map to raw files if skip_trimming set to true in configfile else map to trimmed files

rule map_bowtie2_pe:
    input:
        read1="00.data/00.rawReads/{sample}_1.fastq.gz" if config["skip_trimming"] == "True" else rules.trim_galore_pe.output.fasta_fwd,
        read2="00.data/00.rawReads/{sample}_2.fastq.gz" if config["skip_trimming"] == "True" else rules.trim_galore_pe.output.fasta_rev,
        fasta=expand("00.data/00.Genome/bowtie2_index/{genome}.{id}.bt2",genome=config['genome'], id=["1","2","3","4","rev.1","rev.2"])
    output:
        bam=temp("02.mapping/00.Raw/{sample}.sorted.bam") if config["keep_Raw_BAMs"] != "True" else "02.mapping/00.Raw/{sample}.sorted.bam"
    conda:
        "envs/bowtie2.yaml"
    envmodules:
        "gcc",
        "bowtie2",
        "samtools"
    params:
        basename="02.mapping/00.Raw/{sample}",
        index_path="00.data/00.Genome/bowtie2_index",
        genome="Danio_rerio.GRCz11.dna_sm.primary_assembly",
        additional_params=config["map_bowtie2_pe"]["additional_params"]
    threads: 48
    benchmark: "benchmarks/02.bowtie2_mapping.{sample}.tsv"
    log:
        "logs/02.01.mapping.{sample}.log",
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
        read1="00.data/00.rawReads/{sample}.fastq.gz" if config["skip_trimming"] == "True" else "01.trim/00.trimmedReads/{sample}.SE.fq.gz",
        fasta=expand("00.data/00.Genome/bowtie2_index/{genome}.{id}.bt2",genome=config['genome'], id=["1","2","3","4","rev.1","rev.2"])
    output:
        bam=temp("02.mapping/00.Raw/{sample}.sorted.bam") if config["keep_Raw_BAMs"] != "True" else "02.mapping/00.Raw/{sample}.sorted.bam"
    conda:
        "envs/bowtie2.yaml"
    envmodules:
        "gcc",
        "bowtie2",
        "samtools"
    params:
        basename="02.mapping/00.Raw/{sample}",
        index_path="00.data/00.Genome/bowtie2_index",
        genome="Danio_rerio.GRCz11.dna_sm.primary_assembly",
        additional_params=config["map_bowtie2_se"]["additional_params"]
    threads: 48
    benchmark: "benchmarks/02.bowtie2_mapping.{sample}.tsv"
    log:
        "logs/02.01.mapping.{sample}.log"
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
# (5) #  Remove mitochondrial chromosome                           #
####################################################################

rule remove_mitochondrial_chromosome:
    input:
        bam="02.mapping/00.Raw/{sample}.sorted.bam"
    output:
        intermediate=temp(directory("02.mapping/00.Raw/{sample}")),
        bam_wo_MT=temp("02.mapping/00.Raw/{sample}.no_MT.sorted.bam") if config["keep_BAMs_wo_MT"] != "True" else "02.mapping/00.Raw/{sample}.no_MT.sorted.bam"
    threads: 8
    envmodules:
        "gcc",
        "samtools"
    log:
        "logs/02.03.remove_MT.{sample}.log"
    benchmark: "benchmarks/02.03.remove_MT.{sample}.tsv"
    message:
        "**********Running Rule: remove_mitochondrial_chromosome************"
    shell:
        "scripts/Remove_Mitochondrial_Chromosome.sh {wildcards.sample} {input.bam} {output.intermediate} {output.bam_wo_MT} {threads}"




## TO BE ADDED ##  #### COUNT READS LEFT AFTER REMOVING MITOCHONDRIAL READS <-- done in the "remove_mitochondrial_chromosome" rule

####################################################################
# (6) #  BootStrapping                                             #
####################################################################


#*******#
# (6.0) #   what about test train split?????
#*******#


#*******#
# (6.1) #  Select samples for each round   <-- what about test train split?????
#*******#


rule bootstrap_samples:
    input:
        metadata=SL_array
    output:
        outpath=directory("04.bootstrapping_arrays"),
        expanded_files=expand("04.bootstrapping_arrays/array_{round}.tsv", round=range(1, config["bootstrap_samples"]["rounds"] + 1)),
        #bootstrap_dir=directory("05.bootstrapping"),
        bam=directory(expand("05.bootstrapping/round_{round}/00.mapping", round=range(1, config["bootstrap_samples"]["rounds"] + 1)))
    params:
        rounds=config["bootstrap_samples"]["rounds"],
        samples_per_group=config["bootstrap_samples"]["samples_per_group"],
        seed=config["bootstrap_samples"]["seed"]
    shell:
        """
        python3 scripts/bootstrap.py {input.metadata} {params.rounds} {params.samples_per_group} {output.outpath} --seed {params.seed} {output.bootstrap_dir}
        """


#*******#
# (6.2) #  Downsample samples for all bootstrap roundss
#*******#


rule downsample_bam:
    input:
        raw="02.mapping/00.Raw/{sample}.no_MT.sorted.bam",
        #bam="05.bootstrapping/round_{round}/00.mapping/{sample}.no_MT.sorted.bam"
        bam=rules.remove_mitochondrial_chromosome.output.bam_wo_MT
    output:
        #downsampled_bam="02.mapping/01.DownSampling.{num_reads}/{round}/{sample}.no_MT.sorted.bam"
        downsampled_bam="05.bootstrapping/round_{round}/01.DownSampling.{num_reads}/{sample}.no_MT.sorted.bam"
    conda:
        "envs/bowtie2.yaml"
    envmodules:
        "gcc",
        "samtools"
    params:
        reads=config["DownSample"]["Reads"],
        reads_PE=config["DownSample"]["Reads_PE"],
        reads_SE=config["DownSample"]["Reads_SE"],
        sample_type=lambda wildcards: get_type(wildcards.sample, st)
    threads: config["DownSample"]["threads"]
    benchmark:
        #"benchmarks/02.02.downsample_bam.round_{round}.{sample}.tsv"
        "benchmarks/02.02.downsample_bam.{num_reads}.round_{round}.{sample}.tsv"
    log:
        #"logs/02.02.downsample_bam.{round}.{sample}.log"
        "logs/02.02.downsample_bam.{num_reads}.round_{round}.{sample}.log"
    shell:
        "scripts/downsample_bams.sh {input.bam} {params.reads_PE} {params.reads_SE} {params.sample_type} {output.downsampled_bam} {threads}"

#*******#
# (6.3) #  BAM2BED 
#*******#


rule bamtobed:
    input: 
        #bam="02.mapping/01.DownSampling.{num_reads}/{sample}.no_MT.sorted.bam"
        bam=rules.downsample_bam.output.downsampled_bam
    output:
        #intermediate=temp(directory("03.PeakCalling/00.rawData/DownSampling_{num_reads}/intermediate/{sample}")),
        #bamtobed=temp("03.PeakCalling/00.rawData/DownSampling_{num_reads}/00.bamtobed/{sample}.bamtobed")
        intermediate=temp(directory("05.bootstrapping/round_{round}/02.bamtobed.DS_{num_reads}/intermediate/{sample}")),
        bamtobed=temp("05.bootstrapping/round_{round}/02.bamtobed.DS_{num_reads}/{sample}.bamtobed")
    threads: 32
    params:
        sample_type=lambda wildcards: get_type(wildcards.sample, st)
    envmodules:
        "gcc",
        "bedtools2",
        "samtools"
    log: 
       "logs/03.00.round_{round}.Downsampling_{num_reads}.bamtobed.{sample}.log"
    benchmark: "benchmarks/03.00.round_{round}.Downsampling_{num_reads}.bamtobed.{sample}.tsv"
    message:
        "**********Running Rule: bamtobed************"
    shell:
        "scripts/04.01.GenomeCov_PE_SE.sh {wildcards.sample} {input.bam} {output.intermediate} {output.bamtobed} {threads} {params.sample_type} > {log} 2>&1"



##TO BE ADDED## COUNT READS LEFT AT INTEMEDIATE STEP BEFORE MERGING BOTH PAIRS OF A READ AND RMOVING READS HAVING A LARGE INSERT SIZE

#*******#
# (6.4) #  GenomeCov
#*******#



rule genome_cov:
    input: 
        #bamtobed="03.PeakCalling/00.rawData/DownSampling_{num_reads}/00.bamtobed/{sample}.bamtobed"
        bamtobed=rules.bamtobed.output.bamtobed
    output:
        genomecov_out=temp(
            "05.bootstrapping/round_{round}/03.genomecov.DS_{num_reads}/{sample}.bedgraph"
        ) if config["keep_genomecov_bedgraph"] != "True" else 
            "05.bootstrapping/round_{round}/03.genomecov.DS_{num_reads}/{sample}.bedgraph",
    threads: 8
    log: "logs/03.01.round_{round}.Downsampling_{num_reads}.genome_cov.{sample}.log"
    benchmark: "benchmarks/03.01.round_{round}.Downsampling_{num_reads}.genome_cov.{sample}.tsv"
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
# (6.4* OPTIONAL) #  Filter GenomeCov (Optional)
#*****************#



rule filter_genome_cov:
    input: 
        genomecov=rules.genome_cov.output.genomecov_out
    output:
        genomecov_out_filtered=temp(
            "05.bootstrapping/round_{round}/03.genomecov.DS_{num_reads}/{sample}.filtered_mindepth_{min_depth}.bedgraph"
        ) if config["keep_genomecov_bedgraph"] != "True" else 
            "05.bootstrapping/round_{round}/03.genomecov.DS_{num_reads}/{sample}.filtered_mindepth_{min_depth}.bedgraph"

    threads: 8
    log: "logs/03.01.round_{round}.Downsampling_{num_reads}.genome_cov.filtered_mindepth_{min_depth}.{sample}.log"
    benchmark: "benchmarks/03.01.round_{round}.Downsampling_{num_reads}.genome_cov.filtered_mindepth_{min_depth}.{sample}.tsv"
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



# let the downstream steps use filtered or unfiltered genome_cov file

#*******#
# (6.5) #  PeakCalling
#*******#


seacr_mode=config["peak_calling"]["mode"]
rule peak_calling:
    input:
        bebdgraph_file=rules.filter_genome_cov.output.genomecov_out_filtered 
        if config["peak_calling"]["filter_by_depth"] else 
        rules.genome_cov.output.genomecov_out
        
        #bebdgraph_file="05.bootstrapping/round_{round}/03.genomecov.DS_{num_reads}/{sample}.bedgraph",
    output:
        seacr_out="05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}/{sample}.{seacr_mode}.bed" 
        if config["peak_calling"]["filter_by_depth"] else 
        "05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/SEACR_{seacr_threshold}/{sample}.{seacr_mode}.bed",

    params:
        threshold=config["peak_calling"]["threshold"],
        normalization=config["peak_calling"]["normalization"],
        mode=config["peak_calling"]["mode"],
        output_prefix="05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}/{sample}" 
        if config["peak_calling"]["filter_by_depth"] else 
        "05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/SEACR_{seacr_threshold}/{sample}",
    log:  
       "logs/03.02.round_{round}.Downsampling_{num_reads}.peak_calling.{sample}.filtered_mindepth_{min_depth}.{seacr_mode}.SEACR_{seacr_threshold}.log" 
       if config["peak_calling"]["filter_by_depth"] else 
       "logs/03.02.round_{round}.Downsampling_{num_reads}.peak_calling.{sample}.{seacr_mode}.SEACR_{seacr_threshold}.log"

    shell:
       """
       echo "---- {wildcards.sample}:00:Peak Calling using threshold: {params.threshold} normalization:{params.normalization} mode:{params.mode} ------"
       SEACR_1.3.sh {input.bebdgraph_file} {params.threshold} {params.normalization} {params.mode} {params.output_prefix}  > {log} 2>&1
       """

#### COUNT Peaks Called Per sample
#### COUNT READS In peaks vs elsewhere i.e FRIP


####################################################################
# (7) #  Overlap peaks                                             #
####################################################################


rule overlap_peaks:
    input:
        sample_array="04.bootstrapping_arrays/array_{round}.tsv",
        expanded_files=rules.bootstrap_samples.output.expanded_files,
        peaks=lambda wildcards: expand("05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/SEACR_{seacr_threshold}/{sample}.{seacr_mode}.bed", 
                                        round=wildcards.round, 
                                        seacr_threshold=wildcards.seacr_threshold,
                                        num_reads=wildcards.num_reads, 
                                        sample=get_outfiles(test_samps=None, test_samps_st=None, type="all",array="04.bootstrapping_arrays/array_{}.tsv".format(wildcards.round)),
                                        seacr_mode=config["peak_calling"]["mode"]),
    output:
        merged_regions_binary_data = "05.bootstrapping/05.Overlap_Peaks.DS_{num_reads}/SEACR_{seacr_threshold}.mergepeaks/round_{round}.rawbam.mergedPeaks.binary.tab",
        merged_regions_collapsed = "05.bootstrapping/05.Overlap_Peaks.DS_{num_reads}/SEACR_{seacr_threshold}.mergepeaks/round_{round}.rawbam.mergedPeaks.tab"
    threads: 1
    params:
        #peak_dir="03.PeakCalling/00.rawData/DownSampling_{num_reads}/02.peaks",
        peak_dir="05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/SEACR_{seacr_threshold}",
        bedfile_suffix=expand("{mode}.bed",mode=config["peak_calling"]["mode"]),
        sample_array="04.bootstrapping_arrays/array_{round}.tsv"
    envmodules:
        "gcc",
        "bedtools2"
    log:
        "logs/04.00.round_{round}.DownSampling_{num_reads}.SEACR_{seacr_threshold}.round_{round}.overlap_peaks.log"
    benchmark: "benchmarks/04.00.round_{round}.DownSampling_{num_reads}.SEACR_{seacr_threshold}.round_{round}.overlap_peaks.tsv"
    shell:
        "scripts/MergeSelectedPeakFiles.sh {params.sample_array} {params.peak_dir} {params.bedfile_suffix} {output.merged_regions_binary_data} {output.merged_regions_collapsed}"




use rule overlap_peaks as overlap_peaks_from_filtered_genome_cov with:
    input:
        sample_array="04.bootstrapping_arrays/array_{round}.tsv",
        expanded_files=rules.bootstrap_samples.output.expanded_files,
        peaks=lambda wildcards: expand("05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}/{sample}.{seacr_mode}.bed", 
                                        round=wildcards.round, 
                                        seacr_threshold=wildcards.seacr_threshold,
                                        num_reads=wildcards.num_reads, 
                                        sample=get_outfiles(test_samps=None, test_samps_st=None, type="all",array="04.bootstrapping_arrays/array_{}.tsv".format(wildcards.round)),
                                        seacr_mode=config["peak_calling"]["mode"],
                                        min_depth=wildcards.min_depth),
    output:
        merged_regions_binary_data = "05.bootstrapping/05.Overlap_Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.mergepeaks/round_{round}.rawbam.mergedPeaks.binary.tab",
        merged_regions_collapsed = "05.bootstrapping/05.Overlap_Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.mergepeaks/round_{round}.rawbam.mergedPeaks.tab"
    params:
        peak_dir="05.bootstrapping/round_{round}/04.Peaks.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}",
        bedfile_suffix=expand("{mode}.bed",mode=config["peak_calling"]["mode"]),
        sample_array="04.bootstrapping_arrays/array_{round}.tsv"
    log:
        "logs/04.00.round_{round}.DownSampling_{num_reads}.filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.round_{round}.overlap_peaks.log"
    benchmark: "benchmarks/04.00.round_{round}.DownSampling_{num_reads}.filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.round_{round}.overlap_peaks.tsv"




####################################################################
# (8) #  Generating blacklists                                    #
####################################################################


percentage_threshold=config["generate_suspectList"]["percentage_threshold"]
rule generate_suspectList:
    input:
       raw_peaks= rules.overlap_peaks_from_filtered_genome_cov.output.merged_regions_binary_data if config["peak_calling"]["filter_by_depth"] else rules.overlap_peaks.output.merged_regions_binary_data,
       chr_lengths="00.data/00.Genome/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.fai",
    output:
        SuspectList="05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/round_{round}.SL.bed" 
        if config["peak_calling"]["filter_by_depth"] else 
        "05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/round_{round}.SL.bed",
    threads: 1
    params:
        num_of_samples=config["generate_suspectList"]["num_of_samples"],
        minimum_length=config["generate_suspectList"]["minimum_length"],
        percentage_threshold=config["generate_suspectList"]["percentage_threshold"],
        metadata=SL_array,
        outpath="05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}" 
        if config["peak_calling"]["filter_by_depth"] else 
        "05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}"
    envmodules:
        "gcc",
        "bedtools2",
        "r"
    log:
       "logs/06.Generating_SuspectLists.DS_{num_reads}.filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}.round_{round}.log" 
       if config["peak_calling"]["filter_by_depth"] else 
       "logs/06.Generating_SuspectLists.DS_{num_reads}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}.round_{round}.log"
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



rule Regions_Present_in_all_SL:
    input:
       SL= expand("05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/round_{round}.SL.bed", 
        num_reads=config["DownSample"]["Reads"],
        #round=range(1, config["bootstrap_samples"]["rounds"] + 1),
        round=range(1, 2 + 1),
        seacr_threshold=config["peak_calling"]["threshold"],
        percentage_threshold=config["generate_suspectList"]["percentage_threshold"],
        min_depth=config["peak_calling"]["min_depth"])
        if config["peak_calling"]["filter_by_depth"] else
        expand("05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/round_{round}.SL.bed", 
        num_reads=config["DownSample"]["Reads"],
        #round=range(1, config["bootstrap_samples"]["rounds"] + 1),
        round=range(1, 2 + 1),
        seacr_threshold=config["peak_calling"]["threshold"],
        percentage_threshold=config["generate_suspectList"]["percentage_threshold"]),
       #filt_peaks=rules.overlap_peaks_filtered.output.merged_regions_binary_data,
       chr_lengths="00.data/00.Genome/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.fai",
    output:
        Overlap_SLs="05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/all_rounds.SL.bed" 
        if config["peak_calling"]["filter_by_depth"] else 
        "05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}/all_rounds.SL.bed"
    params:
        input_dir="05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}" 
        if config["peak_calling"]["filter_by_depth"] else 
        "05.bootstrapping/06.Generating_SuspectLists.DS_{num_reads}/SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}"
    threads: 5
    envmodules:
        "gcc",
        "bedtools2",
        "r"
    log:
       "logs/05.00.DownSampled_{num_reads}.filtered_mindepth_{min_depth}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}.log"
       if config["peak_calling"]["filter_by_depth"] else
       "logs/05.00.DownSampled_{num_reads}.SEACR_{seacr_threshold}.SuspectList.prcnt_{percentage_threshold}.log"
    shell:
       """
        scripts/overlap_SL.sh {params.input_dir} {output.Overlap_SLs}  
       """



