rule helpers__fastqc__:
    """Run fastqc over one read file"""
    input:
        fastq="{prefix}.fq.gz",
    output:
        html="{prefix}_fastqc.html",
        zip_="{prefix}_fastqc.zip",
    log:
        "{prefix}_fastqc.log",
    conda:
        "__environment__.yml"
    shell:
        "fastqc --quiet {input} 2> {log} 1>&2"
