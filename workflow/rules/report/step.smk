rule report__step__generate_config_file:
    """Generate the config file for multiqc"""
    output:
        REPORT / "config.yaml",
    log:
        REPORT / "config.log",
    conda:
        "__environment__.yml"
    params:
        chromosome_x=features["sex_chromosomes"][0],
        chromosome_y=features["sex_chromosomes"][1],
    shell:
        """
        echo "samtools_idxstats_xchr: {params.chromosome_x}" >  {output} 2>  {log}
        echo "samtools_idxstats_ychr: {params.chromosome_y}" >> {output} 2>> {log}
        """


rule report__step__reads:
    """Collect all reports for the reads step"""
    input:
        files=rules.reads__fastqc.input,
        config=REPORT / "config.yaml",
    output:
        html=REPORT_STEP / "reads.html",
    log:
        REPORT_STEP / "reads.log",
    conda:
        "__environment__.yml"
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
            {input.files} \
        2> {log} 1>&2
        """


rule report__step__preprocess:
    """Collect all reports for the fastp step"""
    input:
        rules.preprocess__fastp__report.input,
        rules.preprocess__fastp__fastqc.input,
        config=REPORT / "config.yaml",
    output:
        html=REPORT_STEP / "preprocess.html",
    log:
        REPORT_STEP / "preprocess.log",
    conda:
        "__environment__.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title preprocess \
            --force \
            --filename preprocess \
            --outdir {params.dir} \
            --config {input.config} \
            {input} \
        2> {log} 1>&2
        """


rule report__step__quantify:
    """Collect all reports for the star step"""
    input:
        rules.quantify__star__report.input,
        config=REPORT / "config.yaml",
    output:
        html=REPORT_STEP / "quantify.html",
    log:
        REPORT_STEP / "quantify.log",
    conda:
        "__environment__.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title quantify \
            --force \
            --filename quantify \
            --outdir {params.dir} \
            --config {input.config} \
            {input} \
        2> {log} 1>&2
        """


rule report__step:
    """Collect all per step reports for the pipeline"""
    input:
        rules.report__step__reads.output,
        rules.report__step__preprocess.output,
        rules.report__step__quantify.output,
