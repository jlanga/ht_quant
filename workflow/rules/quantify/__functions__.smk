def get_star_out_prefix(wildcards):
    """Get the prefix for STAR output files."""
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    return STAR / f"{sample_id}.{library_id}."
