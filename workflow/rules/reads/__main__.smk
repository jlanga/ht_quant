include: "__functions__.smk"




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
