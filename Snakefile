# path to configuration file for input parameters
configfile: 'config.yaml'
workdir: config["workdir"]
# load required python packages
import pandas as pd

#********************************#
# Define helper python functions #
#********************************#

def get_samples(st):
    """
    return list of samples from samplesheet.tsv 
    """
    return list(st.index)

def get_control(sample_id,st):
    """
    Return the control value for a given sample_id.
    """
    return st.loc[st['sample'] == sample_id, 'control'].values[0]

def get_type(sample_id,st):
    """
    Return the type (Paired End or Single End) for a given sample_id.
    """
    return st.loc[st['sample'] == sample_id, 'type'].values[0]

def get_sampleIDs_with_IgG(st):
    samples_with_IgG=[]
    for sample in st['sample']:
        if pd.notna(st.loc[st['sample'] == sample, 'control'].values[0]):
            samples_with_IgG.append(sample)
    return samples_with_IgG

def IgG_Present(sample,st):
    if 'control' in st.columns:
        if pd.notna(st.loc[st['sample'] == sample, 'control'].values[0]):
            return get_control(sample, st)
    return "none"

def get_trim_infiles(test_samps, test_samps_st, type="SE"): # do we really need this??
    outfiles_SE = []
    outfiles_PE = []
    for sample in test_samps:
        if get_type(sample, test_samps_st) == "PAIRED":
            current_sample=[]
            for id in [1, 2]:
                current_sample.append("00.data/PRJNA738523/{}_{}.fastq.gz".format(sample, id))
            outfiles_PE.append(current_sample)
        elif get_type(sample, test_samps_st) == "SINGLE":
            outfiles_SE.append("00.data/PRJNA738523/{}.fastq.gz".format(sample))
    if type=="SE":  return outfiles_SE
    else: return outfiles_PE

def get_trim_outfiles(test_samps, test_samps_st, type="SE"):
    outfiles_SE = []
    outfiles_PE = []
    for sample in test_samps:
        if get_type(sample, test_samps_st) == "PAIRED":
            for id in [1, 2]:
                outfiles_PE.append("05.Test_SuspectList_with_IgG/01.trimmedReads/{}_{}.PE.fq.gz".format(sample, id))
        elif get_type(sample, test_samps_st) == "SINGLE":
            outfiles_SE.append("05.Test_SuspectList_with_IgG/01.trimmedReads/{}.SE.fq.gz".format(sample))
    if type=="SE":  return outfiles_SE
    else: return outfiles_PE


def get_outfiles(test_samps=None, test_samps_st=None, type="all", array=None):
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


def list_array_files(round, path_prefix="",path_suffix="", downsampled=True):
    if isinstance(round, range):
        rounds = round
    else:
        rounds = [round]
    
    all_files = []
    for rnd in rounds:
        array = f"arrays/array_{rnd}.tsv"
        if downsampled:
            files = [f"{path_prefix}round_{rnd}/{sample}{path_suffix}" for sample in get_outfiles(type="all", array=array)]
            all_files.extend(files)
        else:
            files = [f"{path_prefix}{sample}{path_suffix}" for sample in get_outfiles(type="all", array=array)]
            all_files.extend(files)
    
    return all_files

def check_if_adaptors_provided(st, sample, trim_param):
    if pd.notna(st.loc[st['sample'] == sample, 'Fwd_adaptor'].values[0]):
        a=st.loc[st['sample'] == sample, 'Fwd_adaptor'].values[0]
        trim_param=f"{trim_param} -a '{a}'"
    if  pd.notna(st.loc[st['sample'] == sample, 'Rev_adaptor'].values[0]):
        if get_type(sample,st)=="PAIRED":
            a2=st.loc[st['sample'] == sample, 'Rev_adaptor'].values[0]
            trim_param=f"{trim_param} -a2 '{a2}'"
    return trim_param







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
#ApplySL_array=config["ApplySL_array"]
#test_samps_st = pd.read_table(ApplySL_array).set_index('sample', drop=False)
#test_samps = get_samples(test_samps_st)
#samples_with_IgG=get_sampleIDs_with_IgG(test_samps_st)
#Test_SL_seacr_mode=config["peak_calling"]["mode"]
#Test_SL_seacr_threshold=config["peak_calling"]["threshold"]

## what is this for?
#test_samps_ctrl_st=pd.read_table('05.Test_SuspectList_with_IgG/PRJNA738523.SEACR.array').set_index('sample', drop=False)
## do we need the following??
#SE=get_outfiles(test_samps, test_samps_st, type="SE")
#PE=get_outfiles(test_samps, test_samps_st, type="PE")



SE=get_outfiles(samps, st, type="SE")
PE=get_outfiles(samps, st, type="PE")

##***************************************************##
## Generate SuspectList                              ##
##***************************************************##

# rule Initialize_SuspectPeak_Hunter runs all the other rules required for generating suspectlist
rule Initialize_SuspectPeak_Hunter:
    input:
        #expand("00.data/00.Genome/{genome}.fa.fai" , genome=config["genome"]),
        #expand("01.trim/00.trimmedReads/{sample}_{read}.PE.fq.gz",sample=PE,read=["1","2"]),
        #expand("01.trim/00.trimmedReads/{sample}.SE.fq.gz",sample=SE),
        #expand("01.trim/01.fastqc/{sample}_{id}.PE_fastqc.html",sample=PE,id=["1","2"]),
        #expand("01.trim/01.fastqc/{sample}.SE_fastqc.html",sample=SE),
        expand("02.mapping/00.Raw/{sample}.no_MT.sorted.bam",sample=PE) ,
        expand("02.mapping/00.Raw/{sample}.no_MT.sorted.bam",sample=SE) ,
        expand("02.mapping/01.DownSampling.{num_reads}/{sample}.no_MT.sorted.bam", sample=samps, num_reads=config["DownSample"]["Reads"]),
        expand("03.PeakCalling/00.rawData/DownSampling_{num_reads}/02.peaks/{id}.{seacr_mode}.bed",num_reads=config["DownSample"]["Reads"], id=samps,seacr_mode=config["peak_calling"]["mode"]),
        expand("04.bootstrapping_arrays/array_{round}.tsv", round=list(range(1, config["bootstrap_samples"]["rounds"] + 1))),
        expand("05.Generating_SuspectLists/00.DownSampled_{num_reads}.mergepeaks/round_{round}.rawbam.mergedPeaks.binary.tab",num_reads=config["DownSample"]["Reads"],round=range(1, config["bootstrap_samples"]["rounds"] + 1)),
        #directory(expand("05.Generating_SuspectLists/00.DownSampled_{num_reads}.SuspectList.prcnt_{percentage_threshold}",num_reads=config["DownSample"]["Reads"],percentage_threshold=config["generate_suspectList"]["percentage_threshold"])),
        expand("05.Generating_SuspectLists/00.DownSampled_{num_reads}.SuspectList.prcnt_{percentage_threshold}/round_{round}.SL1.bed",num_reads=config["DownSample"]["Reads"],percentage_threshold=config["generate_suspectList"]["percentage_threshold"],round=list(range(1, config["bootstrap_samples"]["rounds"] + 1))),
        
        
        
        expand("05.Generating_SuspectLists/00.DownSampled_{num_reads}.mergepeaks/{target_group}.rawbam.mergedPeaks.tab",target_group=["TF","active_marks","inactive_marks"],num_reads=config["DownSample"]["Reads"]),
        expand("05.Generating_SuspectLists/00.DownSampled_{num_reads}.SuspectList.prcnt_{percentage_threshold}.Target_Groups/{target_group}.SL1.bed",target_group=["TF","active_marks","inactive_marks"],num_reads=config["DownSample"]["Reads"],percentage_threshold=config["generate_suspectList"]["percentage_threshold"])
        #list_array_files(round=range(1, config["bootstrap_samples"]["rounds"] + 1),path_prefix="03.bootstrapping/01.Downsampled_BAMs.raw/",path_suffix=".sorted.bam"),
        #list_array_files(round=range(1, config["bootstrap_samples"]["rounds"] + 1),path_prefix="03.bootstrapping/01.Downsampled_BAMs.raw/",path_suffix=".no_MT.sorted.bam") 
        #"04.Generating_SuspectLists/00.mergepeaks/rawbam.mergedPeaks.binary.tab",
        #"04.Generating_SuspectLists/00.mergepeaks/filteredbam.mergedPeaks.binary.tab",
        #expand("04.Generating_SuspectLists/SuspectList.prcnt_{percentage_threshold}/BL1_BL2.overlap.out",percentage_threshold=config["generate_suspectList"]["percentage_threshold"]),



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
# (2) #  FastQC and trimming                                       #
####################################################################

rule raw_fastqc:
    input:
        rawread="00.data/00.rawReads/{srr}_{id}.fastq.gz"
    output:
        zip="00.data/01.rawfastqc/{srr}_{id}_fastqc.zip",
        html="00.data/01.rawfastqc/{srr}_{id}_fastqc.html"
    threads:
        32
    params:
        path="00.data/01.rawfastqc/"
    log:
        "logs/00.raw_fastqc_{srr}_{id}.log"
    benchmark: "benchmarks/00.raw_fastqc_{srr}_{id}.tsv"
    shell:
        """
        fastqc {input.rawread} --threads {threads} -o {params.path} 
        """


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

#### COUNT READS LEFT AFTER TRIMMING
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
# (3) #  Prepare Genome for Mapping                                #
####################################################################



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
# (4) #  Mapping                                                   #
####################################################################


rule map_bowtie2_pe:
    input:
        read1=rules.trim_galore_pe.output.fasta_fwd,
        read2=rules.trim_galore_pe.output.fasta_rev,
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
        read1=rules.trim_galore_se.output.fasta,
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


#### COUNT READS LEFT AFTER REMOVING MITOCHONDRIAL READS <-- done in the "remove_mitochondrial_chromosome" rule


####################################################################
# (5) #  BootStrap                                                #
####################################################################


### Select Samples



rule bootstrap_samples:
    input:
        metadata=SL_array
    output:
        outpath=directory("04.bootstrapping_arrays"),
        expanded_files=expand("04.bootstrapping_arrays/array_{round}.tsv", round=list(range(1, config["bootstrap_samples"]["rounds"] + 1)))
    params:
        rounds=config["bootstrap_samples"]["rounds"],
        samples_per_group=config["bootstrap_samples"]["samples_per_group"],
        seed=config["bootstrap_samples"]["seed"]
    shell:
        """
        python3 scripts/bootstrap.py {input.metadata} {params.rounds} {params.samples_per_group} {output.outpath} --seed {params.seed}
        """

### DownSample them

rule downsample_bam:
    input:
        #metadata="03.bootstrapping/00.arrays/array_{round}.tsv",
        #bam="02.mapping/00.Raw/{sample}.no_MT.sorted.bam"
        bam="02.mapping/00.Raw/{sample}.no_MT.sorted.bam"
    output:
        downsampled_bam="02.mapping/01.DownSampling.{num_reads}/{sample}.no_MT.sorted.bam"
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
        "benchmarks/02.02.downsample_bam.{num_reads}.round.{sample}.tsv"
    log:
        #"logs/02.02.downsample_bam.{round}.{sample}.log"
        "logs/02.02.downsample_bam.{num_reads}.{sample}.log"
    shell:
        "scripts/downsample_bams.sh {input.bam} {params.reads_PE} {params.reads_SE} {params.sample_type} {output.downsampled_bam} {threads}"


####################################################################
# (6) #  BAM2BED                                                   #
####################################################################


rule bamtobed:
    input: 
        bam="02.mapping/01.DownSampling.{num_reads}/{sample}.no_MT.sorted.bam"
    output:
        intermediate=temp(directory("03.PeakCalling/00.rawData/DownSampling_{num_reads}/intermediate/{sample}")),
        bamtobed=temp("03.PeakCalling/00.rawData/DownSampling_{num_reads}/00.bamtobed/{sample}.bamtobed")
    threads: 32
    params:
        sample_type=lambda wildcards: get_type(wildcards.sample, st)
    envmodules:
        "gcc",
        "bedtools2",
        "samtools"
    log: 
       "logs/03.00.Downsampling_{num_reads}.bamtobed.{sample}.log"
    benchmark: "benchmarks/03.00.Downsampling_{num_reads}.bamtobed.{sample}.tsv"
    message:
        "**********Running Rule: bamtobed************"
    shell:
        "scripts/04.01.GenomeCov_PE_SE.sh {wildcards.sample} {input.bam} {output.intermediate} {output.bamtobed} {threads} {params.sample_type} > {log} 2>&1"


#### COUNT READS LEFT AT INTEMEDIATE STEP BEFORE MERGING BOTH PAIRS OF A READ AND RMOVING READS HAVING A LARGE INSERT SIZE

####################################################################
# (6) #  GenomeCov                                                 #
####################################################################

rule genome_cov:
    input: 
        bamtobed="03.PeakCalling/00.rawData/DownSampling_{num_reads}/00.bamtobed/{sample}.bamtobed"
    output:
        genomecov_out=temp("03.PeakCalling/00.rawData/DownSampling_{num_reads}/01.genomecov/{sample}.bedgraph") if config["keep_genomecov_bedgraph"] != "True" else "03.PeakCalling/00.rawData/DownSampling_{num_reads}/01.genomecov/{sample}.bedgraph"
    threads: 8
    log: "logs/03.01.Downsampling_{num_reads}.genome_cov.{sample}.log"
    benchmark: "benchmarks/03.01.Downsampling_{num_reads}.genome_cov.{sample}.tsv"
    message:
        "**********Running Rule: genome_cov************"
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


####################################################################
# (7) #  PeakCalling                                               #
####################################################################

seacr_mode=config["peak_calling"]["mode"]
rule peak_calling:
    input:
        bebdgraph_file="03.PeakCalling/00.rawData/DownSampling_{num_reads}/01.genomecov/{id}.bedgraph"
    output:
        seacr_out="03.PeakCalling/00.rawData/DownSampling_{num_reads}/02.peaks/{id}.{seacr_mode}.bed"
    params:
        threshold=config["peak_calling"]["threshold"],
        normalization=config["peak_calling"]["normalization"],
        mode=config["peak_calling"]["mode"],
        output_prefix="03.PeakCalling/00.rawData/DownSampling_{num_reads}/02.peaks/{id}"
    log:  
       "logs/03.02.Downsampling_{num_reads}.peak_calling.{id}.{seacr_mode}.log"
    benchmark: "benchmarks/03.02.Downsampling_{num_reads}.peak_calling.{id}.{seacr_mode}.tsv"
    shell:
       """
       echo "---- {wildcards.id}:00:Peak Calling using threshold: {params.threshold} normalization:{params.normalization} mode:{params.mode} ------"
       SEACR_1.3.sh {input.bebdgraph_file} {params.threshold} {params.normalization} {params.mode} {params.output_prefix}  > {log} 2>&1
       """
#### COUNT Peaks Called Per sample
#### COUNT READS In peaks vs elsewhere i.e FRIP

####################################################################
# (8) #  Overlap peaks                                             #
####################################################################


rule overlap_peaks:
    input:
        peaks=expand("03.PeakCalling/00.rawData/DownSampling_{num_reads}/02.peaks/{sample}.{seacr_mode}.bed", num_reads=config["DownSample"]["Reads"],sample=samps,seacr_mode=config["peak_calling"]["mode"]),
        sample_array="04.bootstrapping_arrays/array_{round}.tsv"
    output:
        merged_regions_binary_data = "05.Generating_SuspectLists/00.DownSampled_{num_reads}.mergepeaks/round_{round}.rawbam.mergedPeaks.binary.tab",
        merged_regions_collapsed = "05.Generating_SuspectLists/00.DownSampled_{num_reads}.mergepeaks/round_{round}.rawbam.mergedPeaks.tab"
    threads: 1
    params:
        peak_dir="03.PeakCalling/00.rawData/DownSampling_{num_reads}/02.peaks",
        bedfile_suffix=expand("{mode}.bed",mode=config["peak_calling"]["mode"])
    envmodules:
        "gcc",
        "bedtools2"
    log:
        "logs/04.00.DownSampling_{num_reads}.round_{round}.overlap_peaks.log"
    benchmark: "benchmarks/04.00.DownSampling_{num_reads}.round_{round}.overlap_peaks.tsv"
    shell:
        "scripts/MergeSelectedPeakFiles.sh {input.sample_array} {params.peak_dir} {params.bedfile_suffix} {output.merged_regions_binary_data} {output.merged_regions_collapsed}"


rule overlap_peaks_by_Target_Groups:
    input:
        peaks=expand("03.PeakCalling/00.rawData/DownSampling_{num_reads}/02.peaks/{sample}.{seacr_mode}.bed", num_reads=config["DownSample"]["Reads"],sample=samps,seacr_mode=config["peak_calling"]["mode"]),
        sample_array="array_{target_group}.tsv"
    output:
        merged_regions_binary_data = "05.Generating_SuspectLists/00.DownSampled_{num_reads}.mergepeaks/{target_group}.rawbam.mergedPeaks.binary.tab",
        merged_regions_collapsed = "05.Generating_SuspectLists/00.DownSampled_{num_reads}.mergepeaks/{target_group}.rawbam.mergedPeaks.tab"
    threads: 1
    params:
        peak_dir="03.PeakCalling/00.rawData/DownSampling_{num_reads}/02.peaks",
        bedfile_suffix=expand("{mode}.bed",mode=config["peak_calling"]["mode"])
    envmodules:
        "gcc",
        "bedtools2"
    log:
        "logs/04.00.DownSampling_{num_reads}.{target_group}.overlap_peaks.log"
    benchmark: "benchmarks/04.00.DownSampling_{num_reads}.{target_group}.overlap_peaks.tsv"
    shell:
        "scripts/MergeSelectedPeakFiles.sh {input.sample_array} {params.peak_dir} {params.bedfile_suffix} {output.merged_regions_binary_data} {output.merged_regions_collapsed}"




####################################################################
# (15) #  Generating blacklists                                    #
####################################################################


percentage_threshold=config["generate_suspectList"]["percentage_threshold"]
rule generate_suspectList_R1:
    input:
       raw_peaks=rules.overlap_peaks.output.merged_regions_binary_data,
       #filt_peaks=rules.overlap_peaks_filtered.output.merged_regions_binary_data,
       chr_lengths="00.data/00.Genome/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.fai",
    output:
        SuspectList="05.Generating_SuspectLists/00.DownSampled_{num_reads}.SuspectList.prcnt_{percentage_threshold}/round_{round}.SL1.bed"
    threads: 1
    params:
        num_of_samples=config["generate_suspectList"]["num_of_samples"],
        minimum_length=config["generate_suspectList"]["minimum_length"],
        percentage_threshold=config["generate_suspectList"]["percentage_threshold"],
        metadata=SL_array,
        outpath="05.Generating_SuspectLists/00.DownSampled_{num_reads}.SuspectList.prcnt_{percentage_threshold}"
    envmodules:
        "gcc",
        "bedtools2",
        "r"
    log:
       "logs/05.00.DownSample_{num_reads}.round_{round}.SuspectList.prcnt_{percentage_threshold}.log"
    message:
       "**********Running Rule: Generating SuspectLists************"
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
       --output_BL_bases={params.outpath}/round_{wildcards.round}.SL1 \
       --ignore_filtered_plots="TRUE" > {log} 2>&1

       sort -k1,1V -k2,2n {output.SuspectList} > {params.outpath}/round_{wildcards.round}.SL1.sorting.bed
       mv {params.outpath}/round_{wildcards.round}.SL1.sorting.bed {output.SuspectList}
       
       """


rule Reagions_Present_in_all_SL:
    input:
       SL=rules.generate_suspectList_R1.output.SuspectList,
       #filt_peaks=rules.overlap_peaks_filtered.output.merged_regions_binary_data,
       chr_lengths="00.data/00.Genome/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.fai",
    output:
        Overlap_SLs="05.Generating_SuspectLists/00.DownSampled_{num_reads}.SuspectList.prcnt_{percentage_threshold}/all_rounds.SL1.bed"
    threads: 1
    params:
        num_of_samples=config["generate_suspectList"]["num_of_samples"],
        minimum_length=config["generate_suspectList"]["minimum_length"],
        percentage_threshold=config["generate_suspectList"]["percentage_threshold"],
        metadata=SL_array,
        outpath="05.Generating_SuspectLists/00.DownSampled_{num_reads}.SuspectList.prcnt_{percentage_threshold}"
    envmodules:
        "gcc",
        "bedtools2",
        "r"
    log:
       "logs/05.00.DownSample_{num_reads}.SuspectList.prcnt_{percentage_threshold}.log"
    message:
       "**********Running Rule: Generating SuspectLists************"
    shell:
       """
        ../../scripts/overlap_SL.sh /work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter/05.Generating_SuspectLists/00.DownSampled_3000000.SuspectList.prcnt_50 /work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter/05.Generating_SuspectLists/00.DownSampled_3000000.SuspectList.prcnt_50/all_rounds.SL1.tsv

       ## the one below is wrogn i guess
       module load gcc bedtools2
       mkdir -p {params.outpath}
       Rscript scripts/Generate_SuspectLists.bootstrap.R \
       --raw_peaks_file={input.raw_peaks} \
       --chr_lengths_file={input.chr_lengths} \
       --metadata={params.metadata} \
       --num_sample={params.num_of_samples} \
       --percentage_threshold={params.percentage_threshold} \
       --min_shared_region_len={params.minimum_length} \
       --output_BL_bases={params.outpath}/round_{wildcards.round}.SL1 \
       --ignore_filtered_plots="TRUE" > {log} 2>&1

       sort -k1,1V -k2,2n {output.SuspectList} > {params.outpath}/round_{wildcards.round}.SL1.sorting.tsv
       mv {params.outpath}/round_{wildcards.round}.SL1.sorting.tsv {output.SuspectList}
       
       """


rule generate_suspectList_for_Target_groups:
    input:
       raw_peaks=rules.overlap_peaks_by_Target_Groups.output.merged_regions_binary_data,
       #filt_peaks=rules.overlap_peaks_filtered.output.merged_regions_binary_data,
       chr_lengths="00.data/00.Genome/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.fai",
    output:
        SuspectList="05.Generating_SuspectLists/00.DownSampled_{num_reads}.SuspectList.prcnt_{percentage_threshold}.Target_Groups/{target_group}.SL1.bed"
    threads: 1
    params:
        num_of_samples=config["generate_suspectList"]["num_of_samples"],
        minimum_length=config["generate_suspectList"]["minimum_length"],
        percentage_threshold=config["generate_suspectList"]["percentage_threshold"],
        metadata=SL_array,
        outpath="05.Generating_SuspectLists/00.DownSampled_{num_reads}.SuspectList.prcnt_{percentage_threshold}.Target_Groups"
    envmodules:
        "gcc",
        "bedtools2",
        "r"
    log:
       "logs/05.00.DownSample_{num_reads}.{target_group}.SuspectList.prcnt_{percentage_threshold}.log"
    message:
       "**********Running Rule: Generating SuspectLists************"
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
       --output_BL_bases={params.outpath}/{wildcards.target_group}.SL1 \
       --ignore_filtered_plots="TRUE" > {log} 2>&1

       sort -k1,1V -k2,2n {output.SuspectList} > {params.outpath}/{wildcards.target_group}.SL1.sorting.bed
       mv {params.outpath}/{wildcards.target_group}.SL1.sorting.bed {output.SuspectList}
       
       """


### Implement Mappability



### Rules on filtered BAMs
####################################################################
# (9) #  Refrerence genome non-repeat region                       #
####################################################################

### in future we would be using regions with high mappability score....... we would keep all reads that do not have atleast 30bp mapped to rehion with good mapability
rule Regions_Complementary_to_Repeats:
    input:
        repeats="00.data/00.Genome/repeatMasker/{genome}.repeats.bed",
        reference="00.data/00.Genome/{genome}.fa.fai"
    output:
        regions_comp_to_repeats=temp("02.mapping/01.Filtered/intermediate/{genome}.Regions_Complementary_to_Repeats.bed"),
    envmodules:
        "gcc",
        "bedtools2"
    log:
       "logs/02.03.{genome}.Regions_Complementary_to_Repeats.log"
    message: 
        "**********Running Rule: Regions_Complementary_to_Repeats************"
    benchmark: "benchmarks/02.03.{genome}.Regions_Complementary_to_Repeats.tsv"
    shell:
       """
       bedtools complement -i {input.repeats} -g <(sort -k1,1 {input.reference}) > {output.regions_comp_to_repeats}
       """


rule GenMap_HighUMapRegions:
    input:
        GenMap="{genome}.bedgraph",
        #GenMap="Danio_rerio.GRCz11.dna_sm.primary_assembly.bedgraph"
        reference_index="00.data/00.Genome/{genome}.fa.fai"
    output:
        LowUMapRegions="{genome}.LowUMapRegions.out",
        HighUMapRgions="{genome}.HighUMapRgions.out",
    params:
        threshold=0.3,
        read_size=150
    envmodules:
        "gcc",
        "bedtools2"
    log:
        "logs/02.03.{genome}.HighUMRegions_genmap.log"
    message:
        "**********Running Rule: GenMap_HighUMapRegions************"
    benchmark: "benchmarks/02.03.{genome}.HighUMRegions_genmap.tsv"
    shell:
       """
       ## run genmap in the same rule
       #awk '{if($4<{params.threshold}) print $1 "\t" $2 "\t" $3+{params.read_size} "\t" $4}' {input.GenMap}  > {output.LowUMapRegions}
       python3 scripts/Extract_LowMapRegions.py {input.GenMap} {params.read_size} {params.threshold} {input.reference_index} > {output.LowUMapRegions}
       bedtools complement -i {sort -k1,1V -k2,2n output.LowUMapRegions} -g <(sort -k1,1V -k2,2n {input.reference_index}) > {output.HighUMapRgions}
       """

####################################################################
# (10) #  Remove multimapped reads and reads overlapping repeats   #
####################################################################

rule Filter_Bams:
    input:
       bam_file="02.mapping/00.Raw/{sample}.sorted.bam",
       regions_comp_to_repeat=expand("02.mapping/01.Filtered/intermediate/{genome}.Regions_Complementary_to_Repeats.bed",genome=genome),
    output:
        intermediate=temp(directory("02.mapping/01.Filtered/intermediate/{sample}")),
        filtered_bam="02.mapping/01.Filtered/{sample}.UM.no_repeats.bam"
    threads: 8
    envmodules:
        "gcc",
        "bedtools2",
        "samtools"
    log: 
       "logs/02.04.Filter_Bams.{sample}.log"
    benchmark: "benchmarks/02.04.Filter_Bams.{sample}.tsv"
    message:
       "**********Running Rule: Filter_Bams************"
    params:
        genome_name=genome,
        genome_path="00.data/00.Genome",
        sample_type=lambda wildcards: get_type(wildcards.sample, st)
    shell:
       """
       scripts/Filter_Bams.sh {wildcards.sample} {params.sample_type} {input.bam_file} {input.regions_comp_to_repeat} {params.genome_path}/{params.genome_name}.fa {output.filtered_bam} {output.intermediate} {threads}
       """


# (10.b) #  Remove multimapped reads and reads overlapping lowMap Regions   #


rule Filter_Bams_GenMap:
    input:
       bam_file="02.mapping/00.Raw/{sample}.sorted.bam",
       regions_comp_to_repeat=expand("00.data/00.Genome/GENMAP/{genome}.HighUMapRegions.{threshold}.out",genome=config["genome"] , threshold=config["GenMap"]["MinMap"]),
    output:
        intermediate=temp(directory("02.mapping/01.Filtered/intermediate/{sample}")),
        filtered_bam="02.mapping/01.Filtered.GenMap/{sample}.UM.no_repeats.bam"
    threads: 8
    envmodules:
        "gcc",
        "bedtools2",
        "samtools"
    log: 
       "logs/02.04.Filter_Bams.GenMap.{sample}.log"
    benchmark: "benchmarks/02.04.Filter_Bams.GenMap.{sample}.tsv"
    message:
       "**********Running Rule: Filter_Bams (GenMap)************"
    params:
        genome_name=config["genome"],
        genome_path="00.data/00.Genome",
        sample_type=lambda wildcards: get_type(wildcards.sample, st)
    shell:
       """
       scripts/Filter_Bams.GenMap.sh {wildcards.sample} {params.sample_type} {input.bam_file} {input.regions_comp_to_repeat} {params.genome_path}/{params.genome_name}.fa {output.filtered_bam} {output.intermediate} {threads}
       """


####################################################################
# (11) #  BAM2BED: Filtered reads                                  #
####################################################################

use rule bamtobed as bamtobed_filtered with:
    input: 
        bam="02.mapping/01.Filtered/{sample}.UM.no_repeats.bam"
    output:
        intermediate=temp(directory("03.PeakCalling/01.Filtered/intermediate.bamtobed/{sample}")),
        bamtobed=temp("03.PeakCalling/01.Filtered/00.bamtobed/{sample}.bamtobed")
    threads: 32
    envmodules:
        "gcc",
        "bedtools2",
        "samtools"
    log: 
       "logs/03.00.filtered.bamtobed.{sample}.log"
    benchmark: "benchmarks/03.00.filtered.bamtobed.{sample}.tsv"
    message:
        "**********Running Rule: bamtobed on filtered Bams************"



####################################################################
# (12) #  GenomeCov: Filtered                                      #
####################################################################

use rule genome_cov as genome_cov_filtered with:
    input: 
        bamtobed="03.PeakCalling/01.Filtered/00.bamtobed/{sample}.bamtobed",
    output:
        genomecov_out=temp("03.PeakCalling/01.Filtered/01.genomecov/{sample}.bedgraph") if config["keep_genomecov_bedgraph"] != "True" else "03.PeakCalling/01.Filtered/01.genomecov/{sample}.bedgraph"
    threads: 8
    log: "logs/03.01.filtered.genome_cov.{sample}.log"
    benchmark: "benchmarks/03.01.filtered.genome_cov.{sample}.tsv"
    message:
        "**********Running Rule: genome_cov on filtered Bams************"


####################################################################
# (13) #  PeakCalling                                               #
####################################################################

use rule peak_calling as peak_calling_filtered with:
    input:
        bebdgraph_file="03.PeakCalling/01.Filtered/01.genomecov/{id}.bedgraph"
    output:
        seacr_out="03.PeakCalling/01.Filtered/02.peaks/{id}.{seacr_mode}.bed"
    params:
        threshold=config["peak_calling"]["threshold"],
        normalization=config["peak_calling"]["normalization"],
        mode=config["peak_calling"]["mode"],
        output_prefix="03.PeakCalling/01.Filtered/02.peaks/{id}"
    log:  
       "logs/03.02.filtered.peak_calling.{id}.{seacr_mode}.log"
    benchmark: "benchmarks/03.02.filtered.peak_calling.{id}.{seacr_mode}.tsv"
    message:
        "**********Running Rule: peak_calling on filtered Bams************"


####################################################################
# (14) #  Overlapping peaks                                        #
####################################################################


use rule overlap_peaks as overlap_peaks_filtered with:
    input:
        peaks=expand("03.PeakCalling/01.Filtered/02.peaks/{sample}.{seacr_mode}.bed", sample=samps,seacr_mode=config["peak_calling"]["mode"])
    output:
        merged_regions_binary_data = "04.Generating_SuspectLists/00.mergepeaks/filteredbam.mergedPeaks.binary.tab",
        merged_regions_collapsed = "04.Generating_SuspectLists/00.mergepeaks/filteredbam.mergedPeaks.tab"
    threads: 1
    log:
        "logs/04.00.filtered.overlap_peaks.log"
    benchmark: "benchmarks/04.00.filtered.overlap_peaks.tsv"
    message:
        "**********Running Rule: overlap_peaks on filtered Bams************"


####################################################################
# (15) #  Generating blacklists                                    #
####################################################################


percentage_threshold=config["generate_suspectList"]["percentage_threshold"]
rule generate_suspectList:
    input:
       #raw_peaks="/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/01_pipeline/06.Peaks.SEACR/threshold_0.001/samples/rep1_targets.stringent.overlap.tab",
       raw_peaks=rules.overlap_peaks.output.merged_regions_binary_data,
       #filt_peaks="/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/01_pipeline/07.Peak_Calling.On_Filtered_Reads/02.SEACR/threshold_0.001/samples/filtered.rep1_targets.stringent.overlap.tab",
       filt_peaks=rules.overlap_peaks_filtered.output.merged_regions_binary_data,
       chr_lengths="00.data/00.Genome/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.fai",
       metadata="00.data/metadata.tsv"
    output:
       outpath=directory("04.Generating_SuspectLists/SuspectList.prcnt_{percentage_threshold}"),
       overlap="04.Generating_SuspectLists/SuspectList.prcnt_{percentage_threshold}/BL1_BL2.overlap.out"
    threads: 1
    params:
        num_of_samples=config["generate_suspectList"]["num_of_samples"],
        minimum_length=config["generate_suspectList"]["minimum_length"],
        percentage_threshold=config["generate_suspectList"]["percentage_threshold"]
    envmodules:
        "gcc",
        "bedtools2",
        "r"
    log:
       "logs/05.00.SuspectList.prcnt_{percentage_threshold}.log"
    message:
       "**********Running Rule: Generating SuspectLists************"
    shell:
       """
       module load gcc bedtools2
       echo "running Generate_SuspectLists.R"
       Rscript scripts/Generate_SuspectLists.R \
       --raw_peaks_file={input.raw_peaks} \
       --filt_peaks_file={input.filt_peaks} \
       --chr_lengths_file={input.chr_lengths} \
       --metadata={input.metadata} \
       --num_sample={params.num_of_samples} \
       --percentage_threshold={params.percentage_threshold} \
       --min_shared_region_len={params.minimum_length} \
       --output_BL_bases={output.outpath}/BL1 \
       --output_BL_bases_filt={output.outpath}/BL2 > {log} 2>&1
       
       echo "sorting and merging BL1 and BL2"
       sort -k 1,1V -k2,2n {output.outpath}/BL1.mns_*.tsv > {output.outpath}/blacklist1.sorted.tsv
       sort -k 1,1V -k2,2n {output.outpath}/BL2.mns_*.tsv > {output.outpath}/blacklist2.sorted.tsv
       bedtools multiinter -i {output.outpath}/blacklist1.sorted.tsv {output.outpath}/blacklist2.sorted.tsv -header -names blacklist1 blacklist2 > {output.overlap}
       
       echo "running Plot_Peaks_Overlap.R"
       Rscript scripts/Plot_Peaks_Overlap.R {output.overlap} {params.minimum_length} {output.outpath}
       """


# collect_stats:
# python3 trimming_summary.py $directory | grep -v "_1.PE*" | awk -F'\t' -v OFS='\t' '{sub(/\.trimming_report/, "", $1); print}' > trimming_summary


#### COUNT READS LEFT AFTER Difefrent Filtering steps

### Get Blacklist <-- ends here.... no ? or should I add apply blacklist part????

