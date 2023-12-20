rule report_library_one:
    """Make a MultiQC report for a single library

    The --dirs option is used to tell MultiQC to append the folder names to
    avoid sample collisions. This is necessary because the same sample name
    can be used for different libraries: reads and fastp.
    """
    input:
        reads=[
            READS / "{sample_id}.{library_id}_1_fastqc.zip",
            READS / "{sample_id}.{library_id}_2_fastqc.zip",
        ],
        fastp=[
            FASTP / "{sample_id}.{library_id}_1_fastqc.zip",
            FASTP / "{sample_id}.{library_id}_2_fastqc.zip",
            FASTP / "{sample_id}.{library_id}_fastp.json",
        ],
        star=get_star_for_library_report,
        config=REPORT / "config.yaml",
    output:
        REPORT_LIBRARY / "{sample_id}.{library_id}.html",
    log:
        REPORT_LIBRARY / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    params:
        library="{sample_id}.{library_id}",
        out_dir=REPORT_LIBRARY,
    shell:
        """
        multiqc \
            --title {params.library} \
            --force \
            --filename {params.library} \
            --outdir {params.out_dir} \
            --dirs \
            --dirs-depth 1 \
            --config {input.config} \
            {input} \
        2> {log} 1>&2
        """


rule report_library_all:
    """Make a MultiQC report for every library"""
    input:
        [
            REPORT_LIBRARY / f"{sample_id}.{library_id}.html"
            for sample_id, library_id in SAMPLE_LIB
        ],


rule report_library:
    """Make all MultiQC reports per library"""
    input:
        rules.report_library_all.input,


localrules:
    report_library_one,
    report_library_all,
    report_library,
