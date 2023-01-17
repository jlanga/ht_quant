# Macro-scale Host Transcriptomics Bioinformatics
This is the repository containing bioinformatic code for macro-scale host transcriptomic data processing. This repository contains all the required information to reproduce the pipeline in any Unix-based environment with a [Conda](https://docs.conda.io/en/latest/) installation. The pipeline is built in [Snakemake](https://snakemake.readthedocs.io/en/stable/), and manages dependencies using conda (or mamba) for reproducibility and deployability.

## Getting started

### What this pipeline requires

This pipeline requires short Illumina paired-end reads (150PE) generated from RNASeq libraries. These are usually produced through poly(A)-enrichment cDNA conversion.

### What this pipeline does

Raw sequencing reads get adapters trimmed and are quality filtered using **Fastp**, before having any leftover ribosomal reads removed using **RiboDetector**. These ribodepleted reads are then mapped to the host genome using **STAR**, which outputs the unmapped reads, sorted BAM file, and count table.

### How this repository is structured

* The 0_Bin directory contains the snakefiles, scripts, and conda environment yml files that are required for running the pipeline.
* The 1_Data directory is were data are processed.
* The 2_Output directory is were the main results are outputed.

## Usage

The easiest way to implement this pipeline is to clone it to your environment and run it within the repository directory.

### 1) Install conda and snakemake

The only prerequisite for running this pipeline is to have conda (or miniconda) and snakemake installed. Miniconda installers for Linux, Mac and Windows operating systems can be found in the following website: https://docs.conda.io/en/latest/miniconda.html

### 2) Clone the repository to the desired directory

Go to the directory in which you want the repository to be cloned, clone the repository using `git clone` and enter into the directory.

``` sh {eval=FALSE}
cd /home/my/desired/directory
git clone https://github.com/3d-omics/Bioinfo_Macro_Host_Transcriptomics.git
cd Bioinfo_Macro_Host_Transcriptomics
```

Note that all paths shown from now on are relative to the path of the repository.

### 3) Place reference genome in the repository

By default the directory in which the reference host genome should be placed is `1_Data/0_Genome/`. This can be modified in the `0_Bin/pipeline.yml` file if you would like the pipeline to use a genome stored elsewhere. Note that the argument `host_genome`in the `yml` file requires the directory where the genome is stored, not the file itself. Hence, the directory can only contain one genome. The required files per genome are the `fasta` file containing the sequence data and `gff` file containing the annotations.

The following code shows an example in which the [chicken genome galgal7](https://www.ncbi.nlm.nih.gov/genome/?term=gallus+gallus) is downloaded from the NCBI.

``` sh {eval=FALSE}
cd 1_Data/0_Genome/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff.gz
cd ../../
```

### 4) Place raw reads in the repository

Raw sequencing reads must be placed in the directory `1_Data/1_Untrimmed` in order for the snakemake pipeline to recognise the input. The data files can be moved to that directory with the following code:

```
mv /home/my/data/*.1.fq.gz 1_Data/1_Untrimmed/
mv /home/my/data/*.2.fq.gz 1_Data/1_Untrimmed/
```

Alternatively, symbolic links can be created to avoid moving or duplicating files:

```
ln -s /home/my/data/*.1.fq.gz 1_Data/1_Untrimmed/
ln -s /home/my/data/*.2.fq.gz 1_Data/1_Untrimmed/
```

### 5) Launch the snakemake pipeline

Once the above-described steps are conducted, the pipeline can be ran using the following script. Note that the `--cluster` argument probably needs to be modified to adjust to your clusters' setup.

```
snakemake \
-s 0_Bin/pipeline.snakefile \
-j 32 \
--cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" \
--use-conda \
--conda-frontend mamba \
--latency-wait 600
```
