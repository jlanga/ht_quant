# Snakemake workflow: `Bioinfo_Macro_Host_Transcriptomics`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/3d-omics/Bioinfo_Macro_Host_Transcriptomics/workflows/Tests/badge.svg?branch=devel)](https://github.com/3d-omics/Bioinfo_Macro_Host_Transcriptomics/actions?query=branch%3Adevel+workflow%3ATests)


A Snakemake workflow for Transcriptome Quantification
##Set up required softwares

the workflow is managed by Snakemake, the reccomended way to install Snakemake is through miniconda
- Install miniconda from terminal (linux)
  - Download the installer script
  ```
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  ```
  - Make the script executable
    ```
    chmod +x Miniconda3-latest-Linux-x86_64.sh
    ```
  - Run the installer script
    ```
    ./Miniconda3-latest-Linux-x86_64.sh
    ```
  - Create a conda environment for Snakemake, you can substitute 'snakemake' with the name you want to assign to the environment
    ```
    conda create -n snakemake
    ```
  - Activate the environment
    ```
    conda activate snakemake
    ```
  - Install snakemake
    ```
    pip install snakemake==7.28.3
    ```
  - Deactivate environment once you are done using snakemake
    ```
    conda deactivate
    ```

## Usage
- Clone the git repository in your terminal
  ```
  git clone git@github.com:3d-omics/Bioinfo_Macro_Host_Transcriptomics.git
  ```
- Change directory to the one you cloned in the previous step
  ```
  cd Bioinfo_Macro_Host_Transcriptomics
  ```
- Test that it works:
  - Make sure you have installed snakemake and run the pipeline with the test
    data, it will download all the necesary software through conda. It should take less than 5 minutes.
  ```
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

  - Run the pipeline
     ```
     snakemake --use-conda --jobs 8 all
     ```
  -  Run the pipeline (slurm users), there is a script called run_slurm in the cloned directory that you can directly use to launch the pipeline on a slurm cluster, you can modify the parameters or direclty execute it as it is
    ```
    ./run_slurm`
    ```


## Features
- FASTQ processing with `fastp`
- Mapping with `STAR`
- SAM/BAM/CRAM processing with `samtools`
- Reports with `multiqc` and `FastQC`
## DAG

![host_transcriptomics_pipeline](./rulegraph.svg?raw=true)

## References

- [`fastp`](https://github.com/OpenGene/fastp)
- [`STAR`](https://github.com/alexdobin/STAR)
- [`samtools`](https://github.com/samtools/samtools)
- [`FastQC`](https://github.com/s-andrews/FastQC)
- [`multiqc`](https://github.com/ewels/MultiQC)
