rule reference__decompress_fasta:
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
        "pigz -dc {input.fa} > {output.fa} 2> {log}"


rule reference__decompress_gtf:
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
        "pigz -dc {input.gtf} > {output.gtf}"


rule reference:
    input:
        rules.reference__decompress_fasta.output,
        rules.reference__decompress_gtf.output,
