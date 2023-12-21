rule __helpers__samtools__index_cram:
    """Generate a cram index"""
    input:
        "{prefix}.cram",
    output:
        "{prefix}.cram.crai",
    log:
        "{prefix}.cram.crai.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools index {input} 2> {log} 1>&2"


rule _helpers__samtools__stats_cram:
    """Compute stats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
        reference=REFERENCE / "genome.fa",
    output:
        tsv="{prefix}.stats.tsv",
    log:
        "{prefix}.stats.log",
    conda:
        "__environment__.yml"
    shell:
        """
        samtools stats \
            --reference {input.reference} \
            {input.cram} \
        > {output.tsv} 2> {log}
        """


rule __helpers__samtools__flagstats_cram:
    """Compute flagstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        txt="{prefix}.flagstats.txt",
    log:
        "{prefix}.flagstats.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools flagstats {input.cram} > {output.txt} 2> {log}"


rule __helpers__samtools__idxstats_cram:
    """Compute idxstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        tsv="{prefix}.idxstats.tsv",
    log:
        "{prefix}.idxstats.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools idxstats {input.cram} > {output.tsv} 2> {log}"
