rule fastp_trim_one:
    """Run fastp on one library"""
    input:
        forward_=READS / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=READS / "{sample_id}.{library_id}_2.fq.gz",
    output:
        forward_=temp(FASTP / "{sample_id}.{library_id}_1.fq.gz"),
        reverse_=temp(FASTP / "{sample_id}.{library_id}_2.fq.gz"),
        unpaired1=temp(FASTP / "{sample_id}.{library_id}_u1.fq.gz"),
        unpaired2=temp(FASTP / "{sample_id}.{library_id}_u2.fq.gz"),
        html=FASTP / "{sample_id}.{library_id}.html",
        json=FASTP / "{sample_id}.{library_id}_fastp.json",
    log:
        FASTP / "{sample_id}.{library_id}.log",
    benchmark:
        FASTP / "{sample_id}.{library_id}.bmk"
    params:
        adapter_forward=get_forward_adapter,
        adapter_reverse=get_reverse_adapter,
        extra=params["fastp"]["extra"],
    threads: 16  # Doesn't work above 16
    resources:
        mem_mb=4 * 1024,
        runtime=24 * 60,
    conda:
        "__environment__.yml"
    shell:
        """
        fastp \
            --in1 {input.forward_} \
            --in2 {input.reverse_} \
            --out1 {output.forward_} \
            --out2 {output.reverse_} \
            --unpaired1 {output.unpaired1} \
            --unpaired2 {output.unpaired2} \
            --html {output.html} \
            --json {output.json} \
            --compression 1 \
            --verbose \
            --trim_poly_g \
            --trim_poly_x \
            --adapter_sequence {params.adapter_forward} \
            --adapter_sequence_r2 {params.adapter_reverse} \
            --thread {threads} \
            {params.extra} \
        2> {log} 1>&2
        """


rule fastp_trim_all:
    """Run fastp over all libraries"""
    input:
        [
            FASTP / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIB
            for end in "1 2 u1 u2".split(" ")
        ],


rule fastp_report_all:
    """Collect fastp reports"""
    input:
        [
            FASTP / f"{sample_id}.{library_id}_fastp.json"
            for sample_id, library_id in SAMPLE_LIB
        ],
        [
            FASTP / f"{sample_id}.{library_id}_{end}_fastqc.zip"
            for sample_id, library_id in SAMPLE_LIB
            for end in "1 2".split(" ")
        ],


rule fastp:
    """Run fastp and collect reports"""
    input:
        rules.fastp_trim_all.input,
        rules.fastp_report_all.input,
