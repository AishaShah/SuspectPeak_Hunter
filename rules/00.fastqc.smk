rule raw_fastqc:
    input:
        rawread="00.data/00.rawReads/{srr}_{id}.fastq.gz"
    output:
        zip="00.data/01.rawfastqc/{srr}_{id}_fastqc.zip",
        html="00.data/01.rawfastqc/{srr}_{id}_fastqc.html"
    threads:
        1
    params:
        path="00.data/01.rawfastqc/"
    log:
        "logs/00.raw_fastqc/00.raw_fastq_{srr}_{id}.log"
    shell:
        """
        fastqc {input.rawread} --threads {threads} -o {params.path} 
        """