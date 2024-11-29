rule multiqc:
    """Collect all reports for the reads step"""
    input:
        reads=[
            READS / f"{sample_id}.{library_id}_{end}_fastqc.zip"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in [1, 2]
        ],
        fastp=[
            FASTP / f"{sample_id}.{library_id}_fastp.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        star=[
            STAR / f"{sample_id}.{library_id}.Log.final.out"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        samtools=[
            STAR / f"{sample_id}.{library_id}.Aligned.sortedByCoord.out.{report}"
            for report in BAM_REPORTS
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        html=RESULTS / "ht_quant.html",
        zip=RESULTS / "ht_quant.zip",
    log:
        RESULTS / "ht_quant.log",
    params:
        extra="--title ht_quant --force",
    wrapper:
        "v5.2.1/bio/multiqc"


rule multiqc__all:
    input:
        RESULTS / "ht_quant.html",
