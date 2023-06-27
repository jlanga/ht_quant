rule fastqc_one:
    """Run fastqc over one read file"""
    input:
        fastq="{prefix}.fq.gz",
    output:
        html="{prefix}_fastqc.html",
        zip_="{prefix}_fastqc.zip",
    log:
        "{prefix}_fastqc.log",
    conda:
        "../envs/fastqc.yml"
    shell:
        "fastqc --quiet {input} 2> {log} 1>&2"
