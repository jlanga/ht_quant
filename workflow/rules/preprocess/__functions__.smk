def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    forward_adapter = samples[
        (samples["sample_id"] == wildcards.sample_id)
        & (samples["library_id"] == wildcards.library_id)
    ]["forward_adapter"].tolist()[0]
    return forward_adapter


def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    reverse_adapter = samples[
        (samples["sample_id"] == wildcards.sample_id)
        & (samples["library_id"] == wildcards.library_id)
    ]["reverse_adapter"].tolist()[0]
    return reverse_adapter
