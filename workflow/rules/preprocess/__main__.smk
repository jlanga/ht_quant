include: "__functions__.smk"
include: "fastp.smk"


rule preprocess:
    input:
        rules.preprocess__fastp.input,
