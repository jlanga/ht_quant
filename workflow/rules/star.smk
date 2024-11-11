rule star__index:
    """Build STAR index of the reference genome"""
    input:
        dna=REFERENCE / (HOST_NAME + ".fa"),
        gtf=REFERENCE / (HOST_NAME + ".gtf"),
    output:
        multiext(
            str(INDEX / HOST_NAME) + "/",
            "chrLength.txt",
            "chrNameLength.txt",
            "chrName.txt",
            "chrStart.txt",
            "exonGeTrInfo.tab",
            "exonInfo.tab",
            "geneInfo.tab",
            "Genome",
            "genomeParameters.txt",
            "Log.out",
            "SA",
            "SAindex",
            "sjdbInfo.txt",
            "sjdbList.fromGTF.out.tab",
            "sjdbList.out.tab",
            "transcriptInfo.tab",
        ),
    log:
        INDEX / (HOST_NAME + ".log"),
    conda:
        "../environments/star.yml"
    cache: True
    params:
        sjdbOverhang=params["quantify"]["star"]["index"]["sjdbOverhang"],
        out_dir=str(INDEX / HOST_NAME) + "/",
    shell:
        """
        STAR \
            --runMode          genomeGenerate \
            --runThreadN       {threads} \
            --genomeFastaFiles {input.dna} \
            --sjdbGTFfile      {input.gtf} \
            --sjdbOverhang     {params.sjdbOverhang} \
            --genomeDir        {params.out_dir} \
        2> {log} 1>&2

        mv Log.out {params.out_dir}
        """


rule star__all:
    input:
        rules.star__index.output,
