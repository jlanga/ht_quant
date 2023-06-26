def get_star_out_prefix(wildcards):
    """Get the prefix for STAR output files."""
    sample = wildcards.sample
    library = wildcards.library
    return STAR / f"{sample}.{library}/{sample}.{library}."
