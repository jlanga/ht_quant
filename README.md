# Macro-scale Host Transcriptomics Bioinformatics
This is the repository containing bioinformatic code for macro-scale host transcriptomic data processing. This repository contains all the required information to reproduce the pipeline in any Unix-based environment with a [Conda](https://docs.conda.io/en/latest/) installation. The pipeline is built on [Snakemake](https://snakemake.readthedocs.io/en/stable/), and manages dependencies using conda (or mamba) for reproducibility and deployability. 

## What this pipeline requires

This pipeline requires short Illumina paired-end reads (150PE) generated from RNASeq libraries. These are usually produced through poly(A)-enrichment cDNA conversion.

## What this pipeline does

Raw sequencing reads get adapters trimmed and are quality filtered using **Fastp**, before having any leftover ribosomal reads removed using **RiboDetector**. These ribodepleted reads are then mapped to the host genome using **STAR**, which outputs the unmapped reads, sorted BAM file, and count table.

## How this repository is structured

* The 0_Bin directory contains the snakefiles, scripts, and conda environment yamls that are required for running the pipeline.
* The 1_Data directory is were data are processed.
* The 2_Results directory is were the main results are outputed.
