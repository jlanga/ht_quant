def get_star_for_library_report(wildcards):
    """Get all star reports for a single library"""
    sample = wildcards.sample
    library = wildcards.library
    files = [STAR / f"{sample}.{library}/{sample}.{library}.Log.final.out"] + [
        STAR / f"{sample}.{library}/{sample}.{library}.{report}"
        for report in BAM_REPORTS
    ]
    return files
