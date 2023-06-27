rule star_index:
    """Build STAR index of the reference genome"""
    input:
        dna=REFERENCE / "genome.fa",
        gtf=REFERENCE / "annotation.gtf",
    output:
        folder=directory("results/star/index"),
    params:
        sjdbOverhang=params["star"]["index"]["sjdbOverhang"],
    conda:
        "../envs/star.yml"
    log:
        REFERENCE / "index.log",
    threads: 24
    resources:
        mem_gb=32,
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


rule star_align_one:
    """Align one library with STAR"""
    input:
        r1=FASTP / "{sample}.{library}_1.fq.gz",
        r2=FASTP / "{sample}.{library}_2.fq.gz",
        index=STAR / "index",
    output:
        bam=temp(
            STAR
            / "{sample}.{library}/{sample}.{library}.Aligned.sortedByCoord.out.bam"
        ),
        counts=STAR / "{sample}.{library}/{sample}.{library}.ReadsPerGene.out.tab",
        out=STAR / "{sample}.{library}/{sample}.{library}.Log.final.out",
    log:
        STAR / "{sample}.{library}/{sample}.{library}.log",
    params:
        out_prefix=get_star_out_prefix,
    conda:
        "../envs/star.yml"
    threads: 24
    resources:
        mem_gb=32,
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


rule star_align_all:
    """Align all libraries with STAR"""
    input:
        [
            STAR / f"{sample}.{library}/{sample}.{library}.ReadsPerGene.out.tab"
            for sample, library in SAMPLE_LIB
        ],


rule star_cram_one:
    """Convert to cram one library

    Note: we use samtools sort when it is already sorted because there is no
    other way to use minimizers on the unmapped fraction.
    """
    input:
        bam=STAR / "{sample}.{library}/{sample}.{library}.Aligned.sortedByCoord.out.bam",
        reference=REFERENCE / "genome.fa",
    output:
        cram=protected(STAR / "{sample}.{library}/{sample}.{library}.cram"),
    log:
        STAR / "{sample}.{library}/{sample}.{library}.cram.log",
    conda:
        "../envs/star.yml"
    threads: 24
    resources:
        mem_gb=32,
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


rule star_cram_all:
    """Convert to cram all the libraries"""
    input:
        [
            STAR / f"{sample}.{library}/{sample}.{library}.cram"
            for sample, library in SAMPLE_LIB
        ],


rule star_create_count_table:
    """Join individual count tables into one"""
    input:
        rules.star_align_all.input,
    output:
        tsv=STAR / "counts.tsv",
    log:
        STAR / "counts.log",
    conda:
        "../envs/star.yml"
    params:
        folder=STAR,
    shell:
        """
        Rscript workflow/scripts/join_star_table.R \
            --input-folder {params.folder} \
            --output-file {output.tsv} \
        2> {log} 1>&2
        """


rule star_report_all:
    """Collect star reports"""
    input:
        [
            STAR / f"{sample}.{library}/{sample}.{library}.Log.final.out"
            for sample, library in SAMPLE_LIB
        ],
        [
            STAR / f"{sample}.{library}/{sample}.{library}.{extension}"
            for sample, library in SAMPLE_LIB
            for extension in BAM_REPORTS
        ],


rule star_all:
    """Run all the star rules for all the libraries"""
    input:
        rules.star_cram_all.input,
        rules.star_report_all.input,
        rules.star_create_count_table.output,


rule star:
    """Run all the star rules"""
    input:
        rules.star_all.input,


localrules:
    star_create_count_table,
