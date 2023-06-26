rule fastqc_one:
    input:
        fastq=READS / "{sample}.{library}_{end}.fq.gz",
    output:
        html=READS / "{sample}.{library}_{end}_fastqc.html",
        zip_=READS / "{sample}.{library}_{end}_fastqc.zip",
    log:
        READS / "{sample}.{library}_{end}_fastqc.log",
    conda:
        "../envs/fastqc.yml"
    shell:
        "fastqc --quiet {input} 2> {log} 1>&2"
