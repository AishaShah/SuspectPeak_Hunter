rule trim_galore_pe:
    input:
        ["00.data/00.rawReads/{sample}_1.fastq.gz", "00.data/00.rawReads/{sample}_2.fastq.gz"],
    output:
        fasta_fwd="01.trim/00.trimmedReads/{sample}_1.fq.gz",
        report_fwd="01.trim/01.Reports/{sample}_1_trimming_report.txt",
        fasta_rev="01.trim/00.trimmedReads/{sample}_2.fq.gz",
        report_rev="01.trim/01.Reports/{sample}_2_trimming_report.txt",
    threads: 1
    params:
        extra="-q 20 --length 20 --stringency 5",
    log:
        "logs/01.trim_galore/{sample}.log",
    wrapper:
        "v3.9.0/bio/trim_galore/pe"
