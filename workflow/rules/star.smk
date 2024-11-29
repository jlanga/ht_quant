# NOTE: do not use wrapper because it uses a directory, and therefore it cannot be cached
rule star__index:
    """Build STAR index of the reference genome"""
    input:
        dna=REFERENCE / f"{HOST_NAME}.fa",
        gtf=REFERENCE / f"{HOST_NAME}.gtf",
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
        INDEX / f"{HOST_NAME}.log",
    conda:
        "../environments/star.yml"
    cache: True
    params:
        sjdbOverhang=params["quantify"]["star"]["index"]["sjdbOverhang"],
        out_dir=lambda w: str(INDEX / HOST_NAME),
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
        """


rule star__index__all:
    input:
        rules.star__index.output,


# NOTE: do not use wrapper because it uses a directory, and therefore it cannot be cached
rule star__align:
    """Align one library with STAR"""
    input:
        forward_=FASTP / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=FASTP / "{sample_id}.{library_id}_2.fq.gz",
        index=multiext(
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
    output:
        bam=STAR / "{sample_id}.{library_id}.Aligned.sortedByCoord.out.bam",
        counts=STAR / "{sample_id}.{library_id}.ReadsPerGene.out.tab",
        out=STAR / "{sample_id}.{library_id}.Log.final.out",
    log:
        STAR / "{sample_id}.{library_id}.log",
    params:
        index_dir=lambda w: str(INDEX / HOST_NAME),
        out_prefix=lambda w: str(STAR / f"{w.sample_id}.{w.library_id}."),
    conda:
        "../environments/star.yml"
    group:
        "quant__{sample_id}.{library_id}"
    shell:
        """
        ulimit -n 90000 2> {log} 1>&2

        STAR \
            --runMode    alignReads \
            --runThreadN {threads} \
            --genomeDir  {params.index_dir} \
            --readFilesIn \
                {input.forward_} \
                {input.reverse_} \
            --outFileNamePrefix {params.out_prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within KeepPairs \
            --readFilesCommand "gzip -cd" \
            --quantMode GeneCounts \
        2>> {log} 1>&2
        """


rule star__align__all:
    input:
        [
            STAR / f"{sample_id}.{library_id}.ReadsPerGene.out.tab"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],


rule star__process__counts:
    input:
        counts=STAR / "{sample_id}.{library_id}.ReadsPerGene.out.tab",
    output:
        STAR / "{sample_id}.{library_id}.tsv",
    log:
        STAR / "{sample_id}.{library_id}.counts.log",
    params:
        tag=lambda w: f"{w.sample_id}.{w.library_id}",
    conda:
        "base"
    shell:
        """
        echo "gene_id\t{params.tag}" > {output} 2> {log}
        ( tail -n+5 {input.counts} | cut -f 1,2 >> {output} ) 2>> {log}
        """


rule star__aggregate:
    """Join individual count tables into one"""
    input:
        [
            STAR / f"{sample_id}.{library_id}.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        tsv=RESULTS / "counts.tsv.gz",
    log:
        RESULTS / "counts.log",
    params:
        subcommand="join",
    wrapper:
        "v5.2.1/utils/csvtk"


rule star__all:
    input:
        rules.star__aggregate.output,
