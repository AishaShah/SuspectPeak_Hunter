rule trimmedReads_fastqc:
    input:
        trimmedread="01.trim/00.trimmedReads/{srr}_{id}.fq.gz"
    output:
        zip="01.trim/01.fastqc/{srr}_{id}_fastqc.zip"
        html="01.trim/01.fastqc/{srr}_{id}_fastqc.html"
    threads:
        8
    params:
        path="01.trim/01.fastqc/"
    log:
        "logs/01.trim_fastqc/00.raw_fastq_{srr}_{id}.log"
    shell:
        """
        fastqc {input.trimmedread} --threads {threads} -o {params.path}
        """
