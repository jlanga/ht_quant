# Snakemake workflow: `Bioinfo_Macro_Host_Transcriptomics`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/3d-omics/Bioinfo_Macro_Host_Transcriptomics/workflows/Tests/badge.svg?branch=devel)](https://github.com/3d-omics/Bioinfo_Macro_Host_Transcriptomics/actions?query=branch%3Adevel+workflow%3ATests)


A Snakemake workflow for Transcriptome Quantification
##Set up required softwares

## Usage
  ```
  #Clone the git repository in your terminal
  git clone git@github.com:3d-omics/Bioinfo_Macro_Host_Transcriptomics.git
  #Change directory to the one you cloned in the previous step
  cd Bioinfo_Macro_Host_Transcriptomics
  #Activate conda environment where you have snakemake
  conda activte Snakemake
  #run the pipeline with the test data, it will download all the necesary software through conda. It should take less than 5 minutes.
  snakemake --use-conda --jobs 8 all
  ```

- Run it with your own data:
  - Edit `config/samples.tsv` and add your samples and where are they located.

  - Run the pipeline
     ```
     snakemake --use-conda --jobs 8 all
     #(slurm users), there is a script called run_slurm in the cloned directory that you can directly use to launch the pipeline on a slurm         cluster, you can modify the parameters or direclty execute it as it is
     ./run_slurm
     ```

## Features

## DAG

![image](https://github.com/3d-omics/Bioinfo_Macro_Host_Transcriptomics/assets/103645443/37f829dd-17e3-4b63-bb42-235714c31520)


## References

