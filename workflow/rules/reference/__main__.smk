rule reference__decompress_fasta__:
    """Decompress the reference genome and save it in the reference folder."""
    input:
        fa=features["dna"],
    output:
        fa=REFERENCE / "genome.fa",
    log:
        REFERENCE / "genome.log",
    conda:
        "__environment__.yml"
    shell:
        "pigz --decompress --stdout {input.fa} > {output.fa} 2> {log}"


rule reference__decompress_gtf__:
    """Decompress the reference annotation and save it in the reference folder."""
    input:
        gtf=features["gtf"],
    output:
        gtf=REFERENCE / "annotation.gtf",
    log:
        REFERENCE / "annotation.log",
    conda:
        "__environment__.yml"
    shell:
        "pigz --decompress --stdout {input.gtf} > {output.gtf}"


rule reference:
    input:
        rules.reference__decompress_fasta__.output,
        rules.reference__decompress_gtf__.output,
