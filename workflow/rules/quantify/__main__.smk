include: "__functions__.smk"
include: "index.smk"
include: "star.smk"


rule quantify:
    input:
        rules.quantify__star.input,
