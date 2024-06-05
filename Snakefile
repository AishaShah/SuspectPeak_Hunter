configfile: 'config.yaml'
import pandas as pd

# map samples to fastqs
def get_samples(st):
    """
    return list of samples from samplesheet.tsv
    """
    return list(st.index)


st = pd.read_table('array.tsv').set_index('sample', drop=False)
samps = get_samples(st)
# Load genome
genome = config["genome"]


rule all:
    input:
        #expand("00.data/01.rawfastqc/{sample}_{read}_fastqc.zip",sample=samps,read=["1","2"]), 
        # expand("01.trim/00.trimmedReads/{sample}_{read}.fq.gz",sample=samps,read=["1","2"]), 
        # expand("01.trim/01.fastqc/{sample}_{read}_fastqc.zip",sample=samps,read=["1","2"]),
        # expand("00.data/00.Genome/bowtie2_index/{genome}.{index}.bt2", genome=genome, index=["1","2","3","4","rev.1","rev.2"]),
        # expand("02.mapping/00.Raw/{sample}.sorted.bam",sample=samps) ,
        # expand("03.PeakCalling/00.rawData/00.bamtobed/{sample}.bamtobed",sample=samps), 
        # expand("03.PeakCalling/00.rawData/01.genomecov/{sample}.bedgraph",sample=samps), 
        # expand("03.PeakCalling/00.rawData/02.peaks/{sample}.{mode}.bed",sample=samps,mode=config["peak_calling"]["mode"]),  ## check this
        #"04.Generating_SuspectLists/00.mergepeaks/rawbam.mergedPeaks.binary.tab" ,
        # expand("02.mapping/01.Filtered/intermediate/{genome}.Regions_Complementary_to_Repeats.bed",genome=genome),
        # expand("02.mapping/01.Filtered/{sample}.UM.no_repeats.SC30.bam",sample=samps) ,
        # expand("03.PeakCalling/01.Filtered/00.bamtobed/{sample}.bamtobed",sample=samps),
        # expand("03.PeakCalling/01.Filtered/01.genomecov/{sample}.bedgraph",sample=samps),
        # expand("03.PeakCalling/01.Filtered/02.peaks/{sample}.{mode}.bed",sample=samps,mode=config["peak_calling"]["mode"]), 
        #expand("04.Generating_SuspectLists/00.mergepeaks/filteredbam.mergedPeaks{extension}",extension=[".tab",".binary.tab"]),
        expand("04.Generating_SuspectLists/SuspectList.prcnt_{percentage_threshold}/BL1_BL2.overlap.out",percentage_threshold=config["generate_suspectList"]["percentage_threshold"]) ## running
        # remove_mitochondrial_chromosome To be tested
#SRR,ID=global_wildcards("00.data/00.rawReads/{srr}_{id}.fastq.gz")




#########################################################################################################


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
    shell:
        """
        fastqc {input.rawread} --threads {threads} -o {params.path} 
        """


#########################################################################################################


rule trim_galore_pe:
    input:
        ["00.data/00.rawReads/{sample}_1.fastq.gz", "00.data/00.rawReads/{sample}_2.fastq.gz"],
    output:
        fasta_fwd="01.trim/00.trimmedReads/{sample}_1.fq.gz",
        report_fwd="01.trim/00.Reports/{sample}_1_trimming_report.txt",
        fasta_rev="01.trim/00.trimmedReads/{sample}_2.fq.gz",
        report_rev="01.trim/00.Reports/{sample}_2_trimming_report.txt",
    threads: 24
    params:
        extra="-q 20 --length 20 --stringency 5",
    log:
        "logs/01.trim_galore_{sample}.log",
    wrapper:
        "v3.9.0/bio/trim_galore/pe"

rule trim_galore_se:
    input:
        "00.data/00.rawReads/{sample}.fastq.gz",
    output:
        fasta="01.trim/00.trimmedReads/{sample}.fq.gz",
        report="01.trim/00.Reports/{sample}_trimming_report.txt",
    params:
        extra="--illumina -q 20 --length 20 --stringency 5",
    log:
        "logs/01.trim_galore_{sample}.log",
    wrapper:
        "v3.10.2/bio/trim_galore/se"

#########################################################################################################


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
    shell:
        """
        fastqc {input.trimmedread} --threads {threads} -o {params.path} 2> {log}
        """


#########################################################################################################


rule genome_bowtie_index:
    input:
        reference="00.data/00.Genome/{genome}.fa"
    output:
        ["00.data/00.Genome/bowtie2_index/{genome}.1.bt2",
        "00.data/00.Genome/bowtie2_index/{genome}.2.bt2",
        "00.data/00.Genome/bowtie2_index/{genome}.3.bt2",
        "00.data/00.Genome/bowtie2_index/{genome}.4.bt2",
        "00.data/00.Genome/bowtie2_index/{genome}.rev.1.bt2",
        "00.data/00.Genome/bowtie2_index/{genome}.rev.2.bt2",
        ]
    params:
        path="00.data/00.Genome/bowtie2_index"
    log:
        "logs/02.00.{genome}_bowtie_index.log",
    threads: 32
    shell:
        """
        bowtie2-build --threads {threads} {input.reference} {params.path}/{genome} 2> {log}
        """


#########################################################################################################


rule map_bowtie2_pe:
    input:
        read1=rules.trim_galore_pe.output.fasta_fwd,
        read2=rules.trim_galore_pe.output.fasta_rev,
        index="00.data/00.Genome/bowtie2_index",
    output:
        bam="02.mapping/00.Raw/{sample}.sorted.bam"
    params:
        basename="02.mapping/00.Raw/{sample}",
        genome="Danio_rerio.GRCz11.dna_sm.primary_assembly"
    threads: 24
    log:
        "logs/02.01.mapping.{sample}.log",
    shell:
       """
       bowtie2 \
       --very-sensitive-local \
       --no-mixed \
       --no-unal \
       --dovetail \
       -X 1000 \
       --threads={threads} \
       -x {input.index}/{params.genome} \
       -1 {input.read1} -2 {input.read2} |\
       samtools sort --threads {threads} -T {params.basename} -o {params.basename}.sorted.bam -
       """ 

rule map_bowtie2_se: ## added
    input:
        read1="01.trim/00.trimmedReads/{sample}.fq.gz",
        index="00.data/00.Genome/bowtie2_index"
    output:
        bam="02.mapping/00.Raw/{sample}.sorted.bam"
    params:
        basename="02.mapping/00.Raw/{sample}",
        genome="Danio_rerio.GRCz11.dna_sm.primary_assembly"
    threads: 24
    log:
        "logs/02.01.mapping.{sample}.log"
    shell:
        """
        bowtie2 --very-sensitive-local --no-unal --threads={threads} -U {input.read1} -x {input.index}/{params.genome} |\
        samtools sort --threads {threads} -T {params.basename} -o {params.basename}.sorted.bam -
        """

rule remove_mitochondrial_chromosome:
    input:
        bam="02.mapping/00.Raw/{sample}.sorted.bam"
    output:
        bam_wo_MT="02.mapping/00.Raw/{sample}.no_MT.sorted.bam",
        intermediate=temp(directory("02.mapping.00.Raw/{sample}"))
    threads: 8
    log:
        "logs/02.02.remove_MT.{sample}.log"
    shell:
        """
        mkdir -p {output.intermediate}
        samtools view --threads {threads} -bS {input.bam} -e 'rname != "MT"' > {output.intermediate}/{wildcards.sample}.no_MT.bam
        samtools sort --threads {threads} {output.intermediate}/{wildcards.sample}.no_MT.bam > {output.bam_wo_MT}
        samtools index -@ {threads} {output.bam_wo_MT}
        """


#########################################################################################################


rule create_genome_index:
    input:
        reference="00.data/00.Genome/{genome}.fa"
    output:
        genome_index="00.data/00.Genome/{genome}.fa.fai"
    threads: 1
    log: 
        "logs/00.Genome_Index.{genome}.log"
    shell:
        """
        samtools faidx {input.reference} {output.genome_index}
        """


#########################################################################################################


rule bamtobed:
    input: 
        bam="02.mapping/00.Raw/{sample}.sorted.bam"
    output:
        intermediate=temp(directory("03.PeakCalling/00.rawData/intermediate/{sample}")),
        bamtobed="03.PeakCalling/00.rawData/00.bamtobed/{sample}.bamtobed"
    threads: 32
    log: 
       "logs/03.00.bamtobed.{sample}.log"
    message:
        "**********Running Rule: bamtobed************"
    shell:
        "scripts/04.01.GenomeCov.sh {wildcards.sample} {input.bam} {output.intermediate} {output.bamtobed} {threads} > {log} 2>&1"


#########################################################################################################


rule genome_cov:
    input: 
        bamtobed="03.PeakCalling/00.rawData/00.bamtobed/{sample}.bamtobed",
    output:
        genomecov_out="03.PeakCalling/00.rawData/01.genomecov/{sample}.bedgraph"
    threads: 8
    log: "logs/03.01.genome_cov.{sample}.log"
    message:
        "**********Running Rule: genome_cov************"
    params:
        genome_name=genome,
        genome_path="00.data/00.Genome"
    shell:
        """
        bedtools genomecov -bg -i {input.bamtobed} -g {params.genome_path}/{params.genome_name}.fa.fai > {output.genomecov_out}
        """
seacr_mode=config["peak_calling"]["mode"]


#########################################################################################################


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
    shell:
       """
       echo "---- {wildcards.id}:00:Peak Calling using threshold: {params.threshold} normalization:{params.normalization} mode:{params.mode} ------"
       SEACR_1.3.sh {input.bebdgraph_file} {params.threshold} {params.normalization} {params.mode} {params.output_prefix}  > {log} 2>&1
       """


#########################################################################################################


rule overlap_peaks:
    input:
        #"03.PeakCalling/00.rawData/02.peaks/{sample}.{seacr_mode}.bed"
        peaks=expand("03.PeakCalling/00.rawData/02.peaks/{sample}.{seacr_mode}.bed", sample=samps,seacr_mode=config["peak_calling"]["mode"])
    output:
        merged_regions_binary_data = "04.Generating_SuspectLists/00.mergepeaks/rawbam.mergedPeaks.binary.tab",
        merged_regions_collapsed = "04.Generating_SuspectLists/00.mergepeaks/rawbam.mergedPeaks.tab"
    threads: 1
    log:
        "logs/04.00.overlap_peaks.log"
    shell:
        """
        bedtools multiinter -i {input.peaks} -header > {output.merged_regions_binary_data}
        scripts/mergePeaks.sh {input.peaks} > {output.merged_regions_collapsed}
        """


#########################################################################################################


### in future we would be using regions with high mappability score....... we would keep all reads that do not have atleast 30bp mapped to rehion with good mapability
rule Regions_Complementary_to_Repeats:
    input:
        repeats="00.data/00.Genome/repeatMasker/{genome}.repeats.bed",
        reference="00.data/00.Genome/{genome}.fa.fai"
    output:
        regions_comp_to_repeats="02.mapping/01.Filtered/intermediate/{genome}.Regions_Complementary_to_Repeats.bed",
    log:
       "logs/02.03.{genome}.Regions_Complementary_to_Repeats.log"
    message: 
        "**********Running Rule: Regions_Complementary_to_Repeats************"
    shell:
       """
       bedtools complement -i {input.repeats} -g <(sort -k1,1 {input.reference}) > {output.regions_comp_to_repeats}
       """


#########################################################################################################


rule Filter_Bams:
    input:
       bam_file="02.mapping/00.Raw/{sample}.sorted.bam",
       regions_comp_to_repeat=expand("02.mapping/01.Filtered/intermediate/{genome}.Regions_Complementary_to_Repeats.bed",genome=genome),
       #regions_comp_to_repeats=rules.Regions_Complementary_to_Repeats.output.regions_comp_to_repeats
    output:
        intermediate=directory("02.mapping/01.Filtered/intermediate/{sample}"),
        filtered_bam="02.mapping/01.Filtered/{sample}.UM.no_repeats.SC30.bam"
    threads: 8
    log: 
       "logs/02.04.Filter_Bams.{sample}.log"
    message:
       "**********Running Rule: Filter_Bams************"
    params:
        genome_name=genome,
        genome_path="00.data/00.Genome"
    shell:
       """
       scripts/03.05.Filter_Bams.sh {wildcards.sample} {input.bam_file} {input.regions_comp_to_repeat} {params.genome_path}/{params.genome_name}.fa {output.filtered_bam} {output.intermediate}
       """

# rule Filter_Bams:
#     input:
#        bam_file="02.mapping/00.Raw/{sample}.sorted.bam",
#        regions_comp_to_repeat=expand("02.mapping/01.Filtered/intermediate/{genome}.Regions_Complementary_to_Repeats.bed",genome=genome),
#     output:
#         intermediate=directory("02.mapping/01.Filtered/intermediate/{sample}"),
#         filtered_bam="02.mapping/01.Filtered/{sample}.UM.no_repeats.SC30.bam"
#     threads: 8
#     log: 
#        "logs/02.04.Filter_Bams.{sample}.log"
#     message:
#        "**********Running Rule: Filter_Bams************"
#     params:
#         genome_name=genome,
#         genome_path="00.data/00.Genome",
#         sample_type=lambda wildcards: get_type(wildcards.sample)
#     shell:
#        """
#        scripts/Filter_Bams.sh {wildcards.sample} {input.bam_file} {input.regions_comp_to_repeat} {params.genome_path}/{params.genome_name}.fa {output.filtered_bam} {output.intermediate} {params.sample_type} {threads}
#        """

#########################################################################################################


use rule bamtobed as bamtobed_filtered with:
    input: 
        bam="02.mapping/01.Filtered/{sample}.UM.no_repeats.SC30.bam"
    output:
        intermediate=temp(directory("03.PeakCalling/01.Filtered/intermediate.bamtobed/{sample}")),
        bamtobed="03.PeakCalling/01.Filtered/00.bamtobed/{sample}.bamtobed"
    threads: 32
    log: 
       "logs/03.00.filtered.bamtobed.{sample}.log"
    message:
        "**********Running Rule: bamtobed on filtered Bams************"



#########################################################################################################

use rule genome_cov as genome_cov_filtered with:
    input: 
        bamtobed="03.PeakCalling/01.Filtered/00.bamtobed/{sample}.bamtobed",
    output:
        genomecov_out="03.PeakCalling/01.Filtered/01.genomecov/{sample}.bedgraph"
    threads: 8
    log: "logs/03.01.filtered.genome_cov.{sample}.log"
    message:
        "**********Running Rule: genome_cov on filtered Bams************"


#########################################################################################################

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
    message:
        "**********Running Rule: peak_calling on filtered Bams************"


#########################################################################################################


use rule overlap_peaks as overlap_peaks_filtered with:
    input:
        peaks=expand("03.PeakCalling/01.Filtered/02.peaks/{sample}.{seacr_mode}.bed", sample=samps,seacr_mode=config["peak_calling"]["mode"])
    output:
        merged_regions_binary_data = "04.Generating_SuspectLists/00.mergepeaks/filteredbam.mergedPeaks.binary.tab",
        merged_regions_collapsed = "04.Generating_SuspectLists/00.mergepeaks/filteredbam.mergedPeaks.tab"
    threads: 1
    log:
        "logs/04.00.filtered.overlap_peaks.log"
    message:
        "**********Running Rule: overlap_peaks on filtered Bams************"


#########################################################################################################


percentage_threshold=config["generate_suspectList"]["percentage_threshold"]
rule generate_suspectList:
    input:
       raw_peaks="/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/01_pipeline/06.Peaks.SEACR/threshold_0.001/samples/rep1_targets.stringent.overlap.tab",
       filt_peaks="/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/01_pipeline/07.Peak_Calling.On_Filtered_Reads/02.SEACR/threshold_0.001/samples/filtered.rep1_targets.stringent.overlap.tab",
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
    log:
       "logs/05.00.SuspectList.prcnt_{percentage_threshold}.log"
    message:
       "**********Running Rule: Generating SuspectLists************"
    shell:
       """
       module load gcc r bedtools2
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
       
       sort -k 1,1V -k2,2n {output.outpath}/BL1.mns_*.tsv > {output.outpath}/blacklist1.sorted.tsv
       sort -k 1,1V -k2,2n {output.outpath}/BL2.mns_*.tsv > {output.outpath}/blacklist2.sorted.tsv
       bedtools multiinter -i {output.outpath}/blacklist1.sorted.tsv {output.outpath}/blacklist2.sorted.tsv -header -names blacklist1 blacklist2 > {output.overlap}
       Rscript scripts/Plot_Peaks_Overlap.R {output.overlap} {params.minimum_length} {output.outpath}
       """



##***************************************************##
## APPLY BLACKLIST  - using IgG                      ##
##***************************************************##


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

test_samps_st = pd.read_table('05.Test_SuspectList_with_IgG/PRJNA738523.samples.array').set_index('sample', drop=False)
#test_samps_st = pd.read_table('05.Test_SuspectList_with_IgG/PRJNA734348.samples.array').set_index('sample', drop=False)

test_samps = get_samples(test_samps_st)
Test_SL_seacr_mode=config["peak_calling"]["mode"]
Test_SL_seacr_threshold=config["peak_calling"]["threshold"]
test_samps_ctrl_st=pd.read_table('05.Test_SuspectList_with_IgG/PRJNA738523.SEACR.array').set_index('sample', drop=False)
samples_with_IgG=get_sampleIDs_with_IgG(test_samps_st)


def get_trim_infiles(test_samps, test_samps_st, type="SE"):
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


SE=get_outfiles(test_samps, test_samps_st, type="SE")
PE=get_outfiles(test_samps, test_samps_st, type="PE")
rule apply_blacklist_with_IgG:
    input:
        expand("05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/Peaks_threshold_{threshold}/testing/{sample}.UM.IgG.{mode}.bed" ,sample=samples_with_IgG,mode=config["peak_calling"]["mode"],threshold=config["peak_calling"]["threshold"])


#### COMMON RULES FOR PEAK CALLING WITH AND WITHOUT SUSPECT LIST REGIONS ######

use rule trim_galore_pe as TestSL_trim_galore_pe with:
    input:
        #["00.data/PRJNA738523/{sample}_1.fastq.gz", "00.data/PRJNA738523/{sample}_2.fastq.gz"],
        ["00.data/PRJNA734348/{sample}_1.fastq.gz", "00.data/PRJNA734348/{sample}_2.fastq.gz"],

    output:
        fasta_fwd="05.Test_SuspectList_with_IgG/01.trimmedReads/{sample}_1.PE.fq.gz",
        report_fwd="05.Test_SuspectList_with_IgG/01.Reports/{sample}_1.PE.trimming_report.txt",
        fasta_rev="05.Test_SuspectList_with_IgG/01.trimmedReads/{sample}_2.PE.fq.gz",
        report_rev="05.Test_SuspectList_with_IgG/01.Reports/{sample}_2.PE.trimming_report.txt",
    log:
        "logs/Test_SuspectList/01.trim_galore_{sample}.log"

use rule trim_galore_se as TestSL_trim_galore_se with:
    input:
        #"00.data/PRJNA738523/{sample}.fastq.gz",
        "00.data/PRJNA734348/{sample}.fastq.gz",
    output:
        fasta="05.Test_SuspectList_with_IgG/01.trimmedReads/{sample}.SE.fq.gz",
        report="05.Test_SuspectList_with_IgG/01.Reports/{sample}.SE.trimming_report.txt"
    log:
        "logs/Test_SuspectList/01.trim_galore_{sample}.log"

use rule map_bowtie2_pe as TestSL_map_bowtie2_pe with:
    input:
        read1=rules.TestSL_trim_galore_pe.output.fasta_fwd,
        read2=rules.TestSL_trim_galore_pe.output.fasta_rev,
        index="00.data/00.Genome/bowtie2_index"
    output:
        bam="05.Test_SuspectList_with_IgG/02.mapping/{sample}.PE.sorted.bam" ### added temp
    params:
        basename="05.Test_SuspectList_with_IgG/02.mapping/{sample}.PE",
        genome="Danio_rerio.GRCz11.dna_sm.primary_assembly"
    threads: 12
    log:
        "logs/Test_SuspectList/02.01.mapping.{sample}.log"


use rule map_bowtie2_se as TestSL_map_bowtie2_se with: ## added
    input:
        read1=rules.TestSL_trim_galore_se.output.fasta,
        index="00.data/00.Genome/bowtie2_index"
    output:
        bam="05.Test_SuspectList_with_IgG/02.mapping/{sample}.SE.sorted.bam"
    params:
        basename="05.Test_SuspectList_with_IgG/02.mapping/{sample}.SE",
        genome="Danio_rerio.GRCz11.dna_sm.primary_assembly"
    threads: 24
    log:
        "logs/Test_SuspectList/02.01.mapping.{sample}.log"

rule TestSL_remove_mitochondrial_chromosome:
    input:
        rules.TestSL_map_bowtie2_pe.output.bam
    output:
        bam_wo_MT="05.Test_SuspectList_with_IgG/02.mapping/{sample}.no_MT.sorted.bam",
        intermediate=temp(directory("05.Test_SuspectList_with_IgG/02.mapping/{sample}"))
    threads: 8
    log:
        "logs/Test_SuspectList/02.02.remove_MT.{sample}.log"
    shell:
        """
        mkdir -p {output.intermediate}
        samtools view --threads {threads} -bS {input} -e 'rname != "MT"' > {output.intermediate}/{wildcards.sample}.no_MT.bam
        samtools sort --threads {threads} {output.intermediate}/{wildcards.sample}.no_MT.bam > {output.bam_wo_MT}
        samtools index -@ {threads} {output.bam_wo_MT}
        """


use rule TestSL_remove_mitochondrial_chromosome as TestSL_remove_mitochondrial_chromosome_SE with:
    input:
        rules.TestSL_map_bowtie2_se.output.bam
    output:
        bam_wo_MT="05.Test_SuspectList_with_IgG/02.mapping/{sample}.no_MT.sorted.bam",
        intermediate=temp(directory("05.Test_SuspectList_with_IgG/02.mapping/{sample}"))
    log:
        "logs/Test_SuspectList/02.02.remove_MT.{sample}.log"


rule Remove_MultiMapped_Reads:
    input:
       bam_file="05.Test_SuspectList_with_IgG/02.mapping/{sample}.no_MT.sorted.bam",
    output:
        intermediate=temp(directory("05.Test_SuspectList_with_IgG/02.mapping/{sample}")),
        filtered_bam="05.Test_SuspectList_with_IgG/02.mapping/{sample}.no_MT.UM.sorted.bam"
    threads: 8
    log: 
        "logs/Test_SuspectList/02.03.remove_MM.{sample}.log"
    message:
       "**********Running Rule: Remove MultiMapped Reads************"
    shell:
       """
       scripts/03.06.Filter_Bams.Remove_MultiMapped.sh {wildcards.sample} {input.bam_file} {output.filtered_bam} {output.intermediate} {threads}
       """


## PEAK CALLING WITHOUT REMOVING SUSPECT LIST REGIONS ######

rule TestSL_bamtobed:
    input: 
        bam="05.Test_SuspectList_with_IgG/02.mapping/{sample}.no_MT.sorted.bam",
        
    output:
        intermediate=temp(directory("05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/intermediate/{sample}")),
        bamtobed="05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/00.bamtobed/{sample}.bamtobed",  ### added temp
    threads: 12
    log: 
       "logs/Test_SuspectList/03.00.wSLR.TestSL_bamtobed.{sample}.log"
    params:
        sample_type=lambda wildcards: get_type(wildcards.sample, test_samps_st)
        #basename="02.mapping/00.Raw/{sample}"
    message:
        "**********Running Rule: bamtobed************"
    shell:
        "scripts/04.01.GenomeCov_PE_SE.sh {wildcards.sample} {input.bam} {output.intermediate} {output.bamtobed} {threads} {params.sample_type} > {log} 2>&1"



use rule TestSL_bamtobed as TestSL_bamtobed_filt with:
    input: 
        bam="05.Test_SuspectList_with_IgG/02.mapping/{sample}.no_MT.UM.sorted.bam",
    output:
        intermediate=temp(directory("05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/intermediate/{sample}.UM")),
        bamtobed="05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/00.bamtobed/{sample}.UM.bamtobed",  ### added temp
    threads: 12
    log: 
       "logs/Test_SuspectList/03.00.wSLR.TestSL_UM.bamtobed.{sample}.log"
    params:
        sample_type=lambda wildcards: get_type(wildcards.sample, test_samps_st)
    message:
        "**********Running Rule: bamtobed_filt************"




use rule genome_cov as TestSL_genomecov with:
    input: 
        bamtobed="05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/00.bamtobed/{sample}{filt}.bamtobed"
    output:
        genomecov_out="05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/01.genomecov/{sample}{filt}.bedgraph" ### added temp
    threads: 12
    log: "logs/Test_SuspectList/03.00.wSLR.TestSL_genomecov.{sample}{filt}.log"
    message:
        "**********Running Rule: genome_cov without Removing SuspectList regions************"
    params:
        genome_name=genome,
        genome_path="00.data/00.Genome"



#######

rule TestSL_Peak_Calling:
    input: 
        bebdgraph_file="05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/01.genomecov/{id}{filt}.bedgraph"
    output:
        seacr_out="05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/Peaks_threshold_{Test_SL_seacr_threshold}/{id}{filt}.{Test_SL_seacr_mode}.bed"
        #seacr_out_IgG="05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/Peaks_threshold_{Test_SL_seacr_threshold}/{id}.IgG.{Test_SL_seacr_mode}.bed"
    threads: 12
    log: 
       "logs/Test_SuspectList/03.01.wSLR.peakcalling.{id}{filt}.{Test_SL_seacr_threshold}.{Test_SL_seacr_mode}.IgG.log"
    params:
        threshold=config["peak_calling"]["threshold"],
        normalization=config["peak_calling"]["normalization"],
        mode=config["peak_calling"]["mode"],
        output_prefix="05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/Peaks_threshold_{Test_SL_seacr_threshold}/{id}{filt}",
    message:
        "**********Running Rule: Peak Calling without Removing SuspectList regions************"
    shell:
       """
       echo "---- {wildcards.id}:00:Peak Calling using threshold: {params.threshold} normalization:{params.normalization} mode:{params.mode} ------"
       SEACR_1.3.sh {input.bebdgraph_file} {params.threshold} {params.normalization} {params.mode} {params.output_prefix}  > {log} 2>&1
       """

# rule TestSL_Peak_Calling_with_IgG:
#     input: 
#         bebdgraph_file="05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/01.genomecov/{id}.{filt}bedgraph"
#     output:
#         seacr_out_IgG="05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/Peaks_threshold_{Test_SL_seacr_threshold}/{id}.{filt}IgG.{Test_SL_seacr_mode}.bed"
#     threads: 12
#     log: 
#        "logs/Test_SuspectList/03.01.wSLR.peakcalling.{id}.{filt}{Test_SL_seacr_threshold}.{Test_SL_seacr_mode}.IgG.log"
#     params:
#         threshold=config["peak_calling"]["threshold"],
#         normalization=config["peak_calling"]["normalization"],
#         mode=config["peak_calling"]["mode"],
#         output_prefix_IgG="05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/Peaks_threshold_{Test_SL_seacr_threshold}/{id}.{filt}IgG",
#         control=lambda wildcards: get_control(wildcards.id, test_samps_st),
#     message:
#         "**********Running Rule: Peak Calling without Removing SuspectList regions************"
#     shell:
#        """
#        echo "---- {wildcards.id}:00:Peak Calling using threshold: {params.control} normalization:norm mode:{params.mode} ------"
#       
#        SEACR_1.3.sh {input.bebdgraph_file} 05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/01.genomecov/{params.control}.bedgraph "norm" {params.mode} {params.output_prefix_IgG}  > {log} 2>&1
#        """
       
rule TestSL_Peak_Calling_with_IgG_testing:
    input: 
        bebdgraph_file="05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/01.genomecov/{id}.UM.bedgraph"
    output:
        seacr_out_IgG="05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/Peaks_threshold_{Test_SL_seacr_threshold}/testing/{id}.UM.IgG.{Test_SL_seacr_mode}.bed"
    threads: 12
    log: 
       "logs/Test_SuspectList/03.01.wSLR.peakcalling.{id}.UM.{Test_SL_seacr_threshold}.{Test_SL_seacr_mode}.IgG.log"
    params:
        threshold=config["peak_calling"]["threshold"],
        normalization=config["peak_calling"]["normalization"],
        mode=config["peak_calling"]["mode"],
        output_prefix_IgG="05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/Peaks_threshold_{Test_SL_seacr_threshold}/testing/{id}.UM.IgG",
        control=lambda wildcards: get_control(wildcards.id, test_samps_st),
    message:
        "**********Running Rule: Peak Calling without Removing SuspectList regions************"
    shell:
       """
       echo "---- {wildcards.id}:00:Peak Calling using threshold: {params.control} normalization:norm mode:{params.mode} ------"
       echo "cmd: {input.bebdgraph_file} 05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/01.genomecov/{params.control}.bedgraph "norm" {params.mode} {params.output_prefix_IgG} "
       SEACR_1.3.sh {input.bebdgraph_file} 05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/01.genomecov/{params.control}.UM.bedgraph "norm" {params.mode} {params.output_prefix_IgG}  > {log} 2>&1
       """
       


# #### PEAK CALLING AFTER REMOVING SUSPECT LIST REGIONS ###### RUN ALL BELOW

rule remove_Suspect_Regions_from_BAMTOBED:
   input: 
      bamtobed="05.Test_SuspectList_with_IgG/03.PeakCalling/01.SL_Regions_not_removed/00.bamtobed/{sample}{filt}.bamtobed",
      suspect_list="04.Generating_SuspectLists/SuspectList.prcnt_{percentage_threshold}/BL1_BL2.overlap.out"
   output: 
      bamtobed_filtered="05.Test_SuspectList_with_IgG/03.PeakCalling/02.SL_Regions_removed.threshold_{percentage_threshold}/00.bamtobed/{sample}{filt}.noSL_Regions.bamtobed" ### added temp
   threads: 12
   log:   "logs/Test_SuspectList/03.00.w0SLR.remove_Suspect_Regions_from_BAMTOBED.{percentage_threshold}.{sample}{filt}.log"
   message:
      "**********Running Rule: Removing reads from bam overlapping Suspect List************"
   shell: 
      """
      echo "---- Removing suspect regions ------"
      bedtools intersect -a {input.bamtobed} -b {input.suspect_list} -v > {output.bamtobed_filtered} 
      """


use rule genome_cov as TestSL_genomecov_SL_Removed with:
    input: 
        bamtobed="05.Test_SuspectList_with_IgG/03.PeakCalling/02.SL_Regions_removed.threshold_{percentage_threshold}/00.bamtobed/{sample}{filt}.noSL_Regions.bamtobed"
    output:
        genomecov_out="05.Test_SuspectList_with_IgG/03.PeakCalling/02.SL_Regions_removed.threshold_{percentage_threshold}/01.genomecov/{sample}{filt}.bedgraph" ### added temp
    threads: 12
    log: "logs/Test_SuspectList/03.00.woSLR.TestSL_genomecov.SL_threshold_{percentage_threshold}.{sample}{filt}.log"
    message:
        "**********Running Rule: genome_cov after removing SuspectList Regions************"
    params:
        genome_name=genome,
        genome_path="00.data/00.Genome"


use rule peak_calling as TestSL_Peak_Calling_SL_Removed with:
    input: 
        bebdgraph_file="05.Test_SuspectList_with_IgG/03.PeakCalling/02.SL_Regions_removed.threshold_{percentage_threshold}/01.genomecov/{id}{filt}.bedgraph"
    output:
        seacr_out="05.Test_SuspectList_with_IgG/03.PeakCalling/02.SL_Regions_removed.threshold_{percentage_threshold}/Peaks_threshold_{Test_SL_seacr_threshold}/{id}{filt}.{Test_SL_seacr_mode}.bed"
    threads: 12
    log: 
       "logs/Test_SuspectList/03.01.woSLR.peakcalling.{percentage_threshold}.{Test_SL_seacr_threshold}.{id}{filt}.{Test_SL_seacr_mode}.log"
    params:
        threshold=config["peak_calling"]["threshold"],
        normalization=config["peak_calling"]["normalization"],
        mode=config["peak_calling"]["mode"],
        output_prefix="05.Test_SuspectList_with_IgG/03.PeakCalling/02.SL_Regions_removed.threshold_{percentage_threshold}/Peaks_threshold_{Test_SL_seacr_threshold}/{id}{filt}"
    message:
        "**********Running Rule: Peak Calling after removing SuspectList Regions************"
