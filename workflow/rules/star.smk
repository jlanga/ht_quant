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
        u1=temp(STAR / "{sample}.{library}/{sample}.{library}.Unmapped.out.mate1"),
        u2=temp(STAR / "{sample}.{library}/{sample}.{library}.Unmapped.out.mate2"),
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
            --outReadsUnmapped Fastx \
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


rule star_import_unmapped_one:
    input:
        u1=STAR / "{sample}.{library}/{sample}.{library}.Unmapped.out.mate1",
        u2=STAR / "{sample}.{library}/{sample}.{library}.Unmapped.out.mate2",
    output:
        bam=temp(STAR / "{sample}.{library}/{sample}.{library}.unmapped.bam"),
    log:
        STAR / "{sample}.{library}/{sample}.{library}.import.log",
    conda:
        "../envs/star.yml"
    threads: 1
    params:
        rg_string="{sample}.{library}",
    resources:
        mem_mb=4 * 1024,
        runtime=24 * 60,
    shell:
        """
        samtools import \
            -1 {input.u1} \
            -2 {input.u2} \
            -R {params.rg_string} \
            -u \
            -o {output.bam} \
        2> {log} 1>&2
        """


rule start_import_unmapped_all:
    input:
        [
            STAR / f"{sample}.{library}/{sample}.{library}.unmapped.bam"
            for sample, library in SAMPLE_LIB
        ],


rule star_join_cram_one:
    """Join mapped and unmapped into one cram file"""
    input:
        mapped_bam=STAR
        / "{sample}.{library}/{sample}.{library}.Aligned.sortedByCoord.out.bam",
        unmapped_bam=STAR / "{sample}.{library}/{sample}.{library}.unmapped.bam",
        reference=REFERENCE / "genome.fa",
    output:
        cram=protected(STAR / "{sample}.{library}/{sample}.{library}.cram"),
    log:
        STAR / "{sample}.{library}/{sample}.{library}.cram.log",
    conda:
        "../envs/star.yml"
    threads: 24
    resources:
        mem_mb=24 * 2048,
        runtime=24 * 60,
    shell:
        """
        (samtools merge \
            -u \
            -o /dev/stdout \
            {input.mapped_bam} \
            {input.unmapped_bam} \
        | samtools sort \
            -l 9 \
            -m 1G \
            -o {output.cram} \
            --output-fmt CRAM \
            --reference {input.reference} \
            -@ {threads} \
            -M \
        ) 2> {log} 1>&2
        """


rule star_cram_all:
    """Convert to cram all the libraries"""
    input:
        [
            STAR / f"{sample}.{library}/{sample}.{library}.cram"
            for sample, library in SAMPLE_LIB
        ],


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


rule star:
    """Run all the star rules"""
    input:
        rules.star_all.input,
