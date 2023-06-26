rule report_generate_config_file:
    output:
        REPORT / "config.yaml",
    log:
        REPORT / "config.log",
    conda:
        "../envs/report.yml"
    params:
        chromosome1=features["sex_chromosomes"][0],
        chromosome2=features["sex_chromosomes"][1],
    shell:
        """
        echo "samtools_idxstats_xchr: {params.chromosome1}" >  {output} 2>  {log}
        echo "samtools_idxstats_ychr: {params.chromosome2}" >> {output} 2>> {log}
        """


rule report_step_reads:
    """Collect all reports for the reads step"""
    input:
        files=rules.reads_fastqc_all.input,
        config=REPORT / "config.yaml",
    output:
        html=REPORT_STEP / "reads.html",
    log:
        REPORT_STEP / "reads.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --filename reads \
            --title reads \
            --force \
            --outdir {params.dir} \
            --config {input.config} \
            {input} \
        2> {log} 1>&2
        """


rule report_step_fastp:
    """Collect all reports for the fastp step"""
    input:
        rules.fastp_report_all.input,
        config=REPORT / "config.yaml",
    output:
        html=REPORT_STEP / "fastp.html",
    log:
        REPORT_STEP / "fastp.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title fastp \
            --force \
            --filename fastp \
            --outdir {params.dir} \
            --config {input.config} \
            {input} \
        2> {log} 1>&2
        """


rule report_step_star:
    """Collect all reports for the star step"""
    input:
        rules.star_report_all.input,
        config=REPORT / "config.yaml",
    output:
        html=REPORT_STEP / "star.html",
    log:
        REPORT_STEP / "star.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title star \
            --force \
            --filename star \
            --outdir {params.dir} \
            --config {input.config} \
            {input} \
        2> {log} 1>&2
        """


rule report_step:
    """Collect all per step reports for the pipeline"""
    input:
        rules.report_step_reads.output,
        rules.report_step_fastp.output,
        rules.report_step_star.output,


localrules:
    report_generate_config_file,
    report_step_reads,
    report_step_fastp,
    report_step_star,
    report_step,
