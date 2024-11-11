include: "preprocess/__functions__.smk"
include: "preprocess/reads.smk"
include: "preprocess/fastp.smk"


rule preprocess__all:
    input:
        rules.preprocess__reads__all.input,
        rules.preprocess__fastp.input,
