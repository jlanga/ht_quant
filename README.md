# Snakemake workflow: `Bioinfo_Macro_Host_Transcriptomics`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/3d-omics/Bioinfo_Macro_Host_Transcriptomics/workflows/Tests/badge.svg?branch=devel)](https://github.com/3d-omics/Bioinfo_Macro_Host_Transcriptomics/actions?query=branch%3Adevel+workflow%3ATests)


##Set up required softwares

## Usage
  ```
  #Clone the git repository in your terminal
  git clone https://github.com/3d-omics/Bioinfo_Macro_Host_Transcriptomics.git
  #Change directory to the one you cloned in the previous step
  cd Bioinfo_Macro_Host_Transcriptomics
  #Activate conda environment where you have snakemake
  conda activte Snakemake
  #run the pipeline with the test data, it will download all the necesary software through conda. It should take less than 5 minutes.
  snakemake --use-conda --jobs 8 all
  ```

- Run it with your own data:
  - Edit `config/samples.tsv` and add your samples and where are they located. Here is an example of the tsv table filled with the information

    ![image](https://github.com/3d-omics/Bioinfo_Macro_Host_Transcriptomics/assets/103645443/bcf67745-9119-498d-a33d-1339ee864246)

  - Edit `config/features.yml` with information regarding the reference you are
    using like in this example.

    ![image](https://github.com/3d-omics/Bioinfo_Macro_Host_Transcriptomics/assets/103645443/195f50ea-eb61-47dd-a650-f91402eca2e3)

  - Edit `config/params.yml` to change the execution of the steps like in this example

    ![image](https://github.com/3d-omics/Bioinfo_Macro_Host_Transcriptomics/assets/103645443/af630e31-c113-4ee6-9408-50870ee54be5)

## Features
- FASTQ processing with `fastp`
- Mapping with `STAR`
- SAM/BAM/CRAM processing with `samtools`
- Reports with `multiqc` and `FastQC`
-
## DAG

![image](https://github.com/3d-omics/Bioinfo_Macro_Host_Transcriptomics/assets/103645443/37f829dd-17e3-4b63-bb42-235714c31520)


## References

- [`fastp`](https://github.com/OpenGene/fastp)
- [`STAR`](https://github.com/alexdobin/STAR)
- [`samtools`](https://github.com/samtools/samtools)
- [`FastQC`](https://github.com/s-andrews/FastQC)
- [`multiqc`](https://github.com/ewels/MultiQC)
