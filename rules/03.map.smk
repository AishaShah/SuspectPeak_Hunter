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
        "logs/03.genome_bowtie_index.log",
    threads: 1
    shell:
        """
        bowtie2-build {input.reference} {params.path}/{genome}
        """
        
    

rule map_bowtie2_pe:
    input:
        read1=rules.trim_galore_pe.output.fasta_fwd,
        read2=rules.trim_galore_pe.output.fasta_rev,
        index="00.data/00.Genome/bowtie2_index",
    output:
        bam="03.mapping/00.Raw/{sample}.sorted.bam"
    params:
        basename="03.mapping/00.Raw/{sample}"
    threads: 16
    log:
        "logs/03.mapping/03.mapping.{sample}.log",
    shell:
       """
       bowtie2 \
       --very-sensitive-local \
       --no-mixed \
       --no-unal \
       --dovetail \
       -X 1000 \
       --threads={threads} \
       -x {input.index} \
       -1 {input.read1} -2 {input.read2} |\
       samtools sort --threads {threads} -T {params.basename} -o {params.basename}.sorted.bam -
       """ 


rule map_boowtie2_se:
    input:
        read=rules.trim_galore_se.output.fasta,
        index="00.data/00.Genome/bowtie2_index",
    output:
        bam="03.mapping/00.Raw/{sample}.sorted.bam"
    params:
        basename="03.mapping/00.Raw/{sample}"
    threads: 16
    log:
        "logs/03.mapping/03.mapping.{sample}.log",
    shell:
       """
       bowtie2 \
       --very-sensitive-local \
       --no-unal \
       --threads={threads} \
       -x {input.index} \
       -U {input.read} |\
       samtools sort --threads {threads} -T {params.basename} -o {params.basename}.sorted.bam -
       """ 