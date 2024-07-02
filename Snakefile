# path to configuration file for input parameters
configfile: 'config.yaml'

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


def get_outfiles(test_samps, test_samps_st, type="SE"):
    outfiles_SE = []
    outfiles_PE = []
    for sample in test_samps:
        if get_type(sample, test_samps_st) == "PAIRED":
                outfiles_PE.append(sample)
        elif get_type(sample, test_samps_st) == "SINGLE":
            outfiles_SE.append(sample)
    if type=="SE":  return list(outfiles_SE)
    else: return list(outfiles_PE)






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
        #"00.data/00.rawReads",
        expand("00.data/00.rawReads/{sample}_{id}.fastq.gz",sample=PE, id=["1","2"]) ,
        expand("00.data/00.rawReads/{sample}.fastq.gz",sample=SE),
        #expand("00.data/00.Genome/{genome}.fa.fai" , genome=config["genome"]),
        #expand("01.trim/00.trimmedReads/{sample}_{read}.PE.fq.gz",sample=PE,read=["1","2"]),
        #expand("01.trim/00.trimmedReads/{sample}.SE.fq.gz",sample=SE),
        #expand("02.mapping/00.Raw/{sample}.sorted.bam",sample=PE) ,
        #expand("02.mapping/00.Raw/{sample}.sorted.bam",sample=SE) ,
        #"04.Generating_SuspectLists/00.mergepeaks/rawbam.mergedPeaks.binary.tab",
        #"04.Generating_SuspectLists/00.mergepeaks/filteredbam.mergedPeaks.binary.tab",
        #expand("04.Generating_SuspectLists/SuspectList.prcnt_{percentage_threshold}/BL1_BL2.overlap.out",percentage_threshold=config["generate_suspectList"]["percentage_threshold"]),



####################################################################
# (1) #  Setup Input Directory <-- create softlinks to inuput data #
####################################################################

rule setup_input_sample_directory:
    input:
        array=SL_array
    output:
        expand("00.data/00.rawReads/{sample}_{id}.fastq.gz",sample=PE, id=["1","2"]) ,
        expand("00.data/00.rawReads/{sample}.fastq.gz",sample=SE) ,
        #directory("01.trim/00.trimmedReads"),
        #directory("02.mapping/00.Raw"),
        directory("00.data/00.rawReads")
    params:
        fq_dir="00.data/00.rawReads",
        fq_trimmed_dir="01.trim/00.trimmedReads",
        bam_dir="02.mapping/00.Raw"
    log:
        "logs/00.setup_input_sample_directory.log"
    shell:
        """
        scripts/Prepare_Input_Sample_directory.sh {input.array} {params.fq_dir} {params.fq_trimmed_dir} {params.bam_dir}
        """

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
    threads: 24
    benchmark: "benchmarks/00.trimming.{sample}.tsv"
    params:
        extra=config["trim_galore_pe"]["additional_params"]
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
    threads: 24
    benchmark: "benchmarks/00.trimming.{sample}.tsv"
    params:
        extra=config["trim_galore_se"]["additional_params"],
    log:
        "logs/01.trim_galore_{sample}.log",
    wrapper:
        "v3.10.2/bio/trim_galore/se"

#### COUNT READS LEFT AFTER TRIMMING

rule trimmedReads_fastqc:
    input:
        trimmedread="01.trim/00.trimmedReads/{srr}_{id}.fq.gz"
    output:
        zip="01.trim/01.fastqc/{srr}_{id}_fastqc.zip",
        html="01.trim/01.fastqc/{srr}_{id}_fastqc.html"
    threads:
        32
    params:
        path="01.trim/01.fastqc/"
    log:
        "logs/01.trimmed_fastq_{srr}_{id}.log"
    benchmark: "benchmarks/01.trimmed_fastq_{srr}_{id}.tsv"
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
        #index="00.data/00.Genome/bowtie2_index",
        fasta=expand("00.data/00.Genome/bowtie2_index/{genome}.{id}.bt2",genome=config['genome'], id=["1","2","3","4","rev.1","rev.2"])
    output:
        bam="02.mapping/00.Raw/{sample}.sorted.bam"
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
    threads: 24
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
        bam="02.mapping/00.Raw/{sample}.sorted.bam"
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
    threads: 24
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
# (5) #  BootStrap                                                #
####################################################################


rule downsample_bam:
    input:
        bam="02.mapping/00.Raw/{sample}.sorted.bam"
    output:
        downsampled_bam="02.mapping/00.Raw/Downsampled/{round}/{sample}.sorted.bam"
    conda:
        "envs/bowtie2.yaml"
    envmodules:
        "gcc",
        "samtools"
    params:
        reads=config["SubSample"]["Reads"],
    threads: 24
    benchmark: "benchmarks/02.02.downsample_bam.{sample}.tsv"
    log:
        "logs/02.02.downsample_bam.{sample}.log"
    shell:
        """
        fraction=$(samtools idxstats {input.bam} | cut -f3 | awk -v ct={params.reads} 'BEGIN {total=0} {total += $1} END {print ct/total}')
        samtools view -b -s ${fraction} {input.bam} > {output.downsampled_bam}
        """


####################################################################
# (5) #  Remove mitochondrial chromosome                           #
####################################################################

rule remove_mitochondrial_chromosome:
    input:
        bam="02.mapping/00.Raw/{sample}.sorted.bam"
    output:
        intermediate=temp(directory("02.mapping/00.Raw/{sample}")),
        bam_wo_MT=temp("02.mapping/00.Raw/{sample}.no_MT.sorted.bam") if config["keep_filtered_BAMs"] != "True" else "02.mapping/00.Raw/{sample}.no_MT.sorted.bam"
    threads: 8
    envmodules:
        "gcc",
        "samtools"
    log:
        "logs/02.03.remove_MT.{sample}.log"
    benchmark: "benchmarks/02.03.remove_MT.{sample}.tsv"
    shell:
        """
        mkdir -p {output.intermediate}
        samtools view --threads {threads} -bS {input.bam} -e 'rname != "MT"' > {output.intermediate}/{wildcards.sample}.no_MT.bam
        samtools sort --threads {threads} {output.intermediate}/{wildcards.sample}.no_MT.bam > {output.bam_wo_MT}
        samtools index -@ {threads} {output.bam_wo_MT}
        """

#### COUNT READS LEFT AFTER REMOVING MITOCHONDRIAL READS

####################################################################
# (6) #  BAM2BED                                                   #
####################################################################


rule bamtobed:
    input: 
        bam="02.mapping/00.Raw/{sample}.no_MT.sorted.bam"
    output:
        intermediate=temp(directory("03.PeakCalling/00.rawData/intermediate/{sample}")),
        bamtobed=temp("03.PeakCalling/00.rawData/00.bamtobed/{sample}.bamtobed")
    threads: 32
    params:
        sample_type=lambda wildcards: get_type(wildcards.sample, st)
    envmodules:
        "gcc",
        "bedtools2",
        "samtools"
    log: 
       "logs/03.00.bamtobed.{sample}.log"
    benchmark: "benchmarks/03.00.bamtobed.{sample}.tsv"
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
        bamtobed="03.PeakCalling/00.rawData/00.bamtobed/{sample}.bamtobed"
    output:
        genomecov_out=temp("03.PeakCalling/00.rawData/01.genomecov/{sample}.bedgraph") if config["keep_genomecov_bedgraph"] != "True" else "03.PeakCalling/00.rawData/01.genomecov/{sample}.bedgraph"
    threads: 8
    log: "logs/03.01.genome_cov.{sample}.log"
    benchmark: "benchmarks/03.01.genome_cov.{sample}.tsv"
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
        bebdgraph_file="03.PeakCalling/00.rawData/01.genomecov/{id}.bedgraph"
    output:
        seacr_out="03.PeakCalling/00.rawData/02.peaks/{id}.{seacr_mode}.bed"
    params:
        threshold=config["peak_calling"]["threshold"],
        normalization=config["peak_calling"]["normalization"],
        mode=config["peak_calling"]["mode"],
        output_prefix="03.PeakCalling/00.rawData/02.peaks/{id}"
    log:  
       "logs/03.02.peak_calling.{id}.{seacr_mode}.log"
    benchmark: "benchmarks/03.02.peak_calling.{id}.{seacr_mode}.tsv"
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
        peaks=expand("03.PeakCalling/00.rawData/02.peaks/{sample}.{seacr_mode}.bed", sample=samps,seacr_mode=config["peak_calling"]["mode"])
    output:
        merged_regions_binary_data = "04.Generating_SuspectLists/00.mergepeaks/rawbam.mergedPeaks.binary.tab",
        merged_regions_collapsed = "04.Generating_SuspectLists/00.mergepeaks/rawbam.mergedPeaks.tab"
    threads: 1
    envmodules:
        "gcc",
        "bedtools2"
    log:
        "logs/04.00.overlap_peaks.log"
    benchmark: "benchmarks/04.00.overlap_peaks.tsv"
    shell:
        """
        bedtools multiinter -i {input.peaks} -header > {output.merged_regions_binary_data}
        scripts/mergePeaks.sh {input.peaks} > {output.merged_regions_collapsed}
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

