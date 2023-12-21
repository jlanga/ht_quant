def get_star_for_library_report(wildcards):
    """Get all star reports for a single library"""
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    files = [STAR / f"{sample_id}.{library_id}.Log.final.out"] + [
        STAR / f"{sample_id}.{library_id}.{report}" for report in BAM_REPORTS
    ]
    return files
