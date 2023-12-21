rule _quantify__star__index:
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


rule _quantify__star__align:
    """Align one library with STAR"""
    input:
        r1=FASTP / "{sample_id}.{library_id}_1.fq.gz",
        r2=FASTP / "{sample_id}.{library_id}_2.fq.gz",
        index=INDEX,
    output:
        bam=temp(STAR / "{sample_id}.{library_id}.Aligned.sortedByCoord.out.bam"),
        counts=STAR / "{sample_id}.{library_id}.ReadsPerGene.out.tab",
        out=STAR / "{sample_id}.{library_id}.Log.final.out",
    log:
        STAR / "{sample_id}.{library_id}.log",
    params:
        out_prefix=get_star_out_prefix,
    conda:
        "__environment__.yml"
    threads: 24
    resources:
        mem_mb=32 * 1024,
        runtime=24 * 60,
    shell:
        """
        ulimit -n 90000 2> {log} 1>&2

        STAR \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {input.index} \
            --readFilesIn {input.r1} {input.r2} \
            --outFileNamePrefix {params.out_prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within KeepPairs \
            --readFilesCommand "gzip -cd" \
            --quantMode GeneCounts \
        2>> {log} 1>&2
        """


rule quantify__star__align:
    """Align all libraries with STAR"""
    input:
        [
            STAR / f"{sample_id}.{library_id}.ReadsPerGene.out.tab"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],


rule _quantify__star__bam_to_cram:
    """Convert to cram one library

    Note: we use samtools sort when it is already sorted because there is no
    other way to use minimizers on the unmapped fraction.
    """
    input:
        bam=STAR / "{sample_id}.{library_id}.Aligned.sortedByCoord.out.bam",
        reference=REFERENCE / "genome.fa",
    output:
        cram=STAR / "{sample_id}.{library_id}.cram",
    log:
        STAR / "{sample_id}.{library_id}.cram.log",
    conda:
        "__environment__.yml"
    threads: 24
    resources:
        mem_mb=32 * 1024,
        runtime=24 * 60,
    shell:
        """
        samtools sort \
            -l 9 \
            -m 1G \
            -o {output.cram} \
            --output-fmt CRAM \
            --reference {input.reference} \
            -@ {threads} \
            -M \
            {input.bam} \
        2> {log} 1>&2
        """


rule quantify__star__bam_to_cram:
    """Convert to cram all the libraries"""
    input:
        [
            STAR / f"{sample_id}.{library_id}.cram"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],


rule _quantify__star__aggregate_counts:
    """Join individual count tables into one"""
    input:
        rules.quantify__star__align.input,
    output:
        tsv=QUANTIFY / "counts.tsv",
    log:
        QUANTIFY / "counts.log",
    conda:
        "__environment__.yml"
    params:
        folder=STAR,
    shell:
        """
        Rscript workflow/scripts/join_star_table.R \
            --input-folder {params.folder} \
            --output-file {output.tsv} \
        2> {log} 1>&2
        """


rule quantify__star__report:
    """Get all reports for star"""
    input:
        logs=[
            STAR / f"{sample_id}.{library_id}.Log.final.out"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
        + [
            STAR / f"{sample_id}.{library_id}.{report}"
            for report in BAM_REPORTS
            for sample_id, library_id in SAMPLE_LIBRARY
        ],


rule quantify__star:
    """Run all the star rules"""
    input:
        rules.quantify__star__bam_to_cram.input,
        rules._quantify__star__aggregate_counts.output,
