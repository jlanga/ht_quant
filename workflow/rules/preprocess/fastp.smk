rule __preprocess__fastp__trim:
    """Run fastp on one library"""
    input:
        forward_=READS / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=READS / "{sample_id}.{library_id}_2.fq.gz",
    output:
        forward_=temp(FASTP / "{sample_id}.{library_id}_1.fq.gz"),
        reverse_=temp(FASTP / "{sample_id}.{library_id}_2.fq.gz"),
        unpaired1=temp(FASTP / "{sample_id}.{library_id}_u1.fq.gz"),
        unpaired2=temp(FASTP / "{sample_id}.{library_id}_u2.fq.gz"),
        html=FASTP / "{sample_id}.{library_id}_fastp.html",
        json=FASTP / "{sample_id}.{library_id}_fastp.json",
    log:
        FASTP / "{sample_id}.{library_id}.log",
    benchmark:
        FASTP / "{sample_id}.{library_id}.bmk"
    params:
        adapter_forward=get_forward_adapter,
        adapter_reverse=get_reverse_adapter,
    threads: 24
    resources:
        mem_mb=4 * 1024,
        runtime=1 * 60,
    conda:
        "__environment__.yml"
    shell:
        """
        fastp \
            --in1 {input.forward_} \
            --in2 {input.reverse_} \
            --out1 >(pigz --fast > {output.forward_}) \
            --out2 >(pigz --fast > {output.reverse_}) \
            --unpaired1 >(pigz --fast > {output.unpaired1}) \
            --unpaired2 >(pigz --fast > {output.unpaired2}) \
            --html {output.html} \
            --json {output.json} \
            --compression 1 \
            --verbose \
            --trim_poly_g \
            --trim_poly_x \
            --adapter_sequence {params.adapter_forward} \
            --adapter_sequence_r2 {params.adapter_reverse} \
            --thread {threads} \
        2> {log} 1>&2
        """


rule preprocess__fastp__trim:
    """Run fastp over all libraries"""
    input:
        [
            FASTP / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in "1 2 u1 u2".split(" ")
        ],


rule preprocess__fastp__fastqc:
    """Run fastqc over all libraries"""
    input:
        [
            FASTP / f"{sample_id}.{library_id}_{end}_fastqc.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in ["1", "2"]
            for extension in ["html", "zip"]
        ],


rule preprocess__fastp__report:
    """Collect fastp and fastqc reports"""
    input:
        [
            FASTP / f"{sample_id}.{library_id}_fastp.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        rules.preprocess__fastp__fastqc.input,


rule preprocess__fastp:
    """Run fastp and collect reports"""
    input:
        rules.preprocess__fastp__trim.input,
        rules.preprocess__fastp__report.input,
