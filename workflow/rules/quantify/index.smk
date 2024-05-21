rule quantify__index__:
    """Build STAR index of the reference genome"""
    input:
        dna=REFERENCE / "genome.fa",
        gtf=REFERENCE / "annotation.gtf",
    output:
        folder=directory(INDEX),
    params:
        sjdbOverhang=params["quantify"]["star"]["index"]["sjdbOverhang"],
    conda:
        "__environment__.yml"
    log:
        REFERENCE / "index.log",
    threads: 24
    resources:
        mem_mb=32 * 1024,
        runtime=24 * 60,
    shell:
        """
        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output.folder} \
            --genomeFastaFiles {input.dna} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.sjdbOverhang} \
        2> {log} 1>&2
        """

rule quantify__index:
    input:
        directory(INDEX)
