include: "__functions__.smk"


rule reads__link__:
    """Make a link to the original file, with a prettier name than default"""
    input:
        forward_=get_forward,
        reverse_=get_reverse,
    output:
        forward_=READS / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=READS / "{sample_id}.{library_id}_2.fq.gz",
    log:
        READS / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.forward_}) {output.forward_}
        ln --symbolic $(readlink --canonicalize {input.reverse_}) {output.reverse_}
        """


rule reads__link:
    input:
        [
            READS / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],


rule reads__fastqc:
    """Run fastqc on all raw reads"""
    input:
        [
            READS / f"{sample_id}.{library_id}_{end}_fastqc.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in ["1", "2"]
            for extension in ["html", "zip"]
        ],


rule reads:
    """Link all reads and run fastqc on them"""
    input:
        rules.reads__link.input,
        rules.reads__fastqc.input,
