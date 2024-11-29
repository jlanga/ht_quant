def compose_adapters(wildcards):
    """Compose the forward and reverse adapter line for fastp"""
    forward_, reverse_ = (
        samples[
            (samples["sample_id"] == wildcards.sample_id)
            & (samples["library_id"] == wildcards.library_id)
        ][["forward_adapter", "reverse_adapter"]]
        .iloc[0]
        .tolist()
    )

    return f"--adapter_sequence {forward_} --adapter_sequence_r2 {reverse_}"
