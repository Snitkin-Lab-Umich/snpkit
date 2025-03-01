[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10](https://img.shields.io/badge/python-3.1-blue.svg)](https://www.python.org/downloads/release/python-310/)

# SNPKIT - Microbial Variant Calling and Diagnostics toolkit.

## Fork Notice:

This repository is the actively maintained version of the snpkit pipeline, originally forked from [https://github.com/alipirani88/snpkit.wiki.git](https://github.com/alipirani88/snpkit/) on 12/30/2024 by Kyle Gontjes. Please use this repository instead of the original.

This version will be maintained until the Snakemake implementation of snpkit (https://github.com/Snitkin-Lab-Umich/snpkit-smk) is complete.  

## Description

SNPKIT is a variant detection workflow that can be easily deployed for infectious disease outbreak investigations and other clinical microbiology projects. The workflow takes Illumina fastq reads and an annotated reference genome as input, calls variants using SAMTOOLS, GATK, and Freebayes, and generates an annotated SNP/Indel matrix for variant diagnostics and visualizations. 

While this pipeline has downstream scripts to generate a Gubbins recombination-filtered phylogenetic tree using the identified variants, we suggest that phylogenetic trees are reconstructed using our updated phylogenetics pipeline: https://github.com/Snitkin-Lab-Umich/phylokit.   

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Input](#input)

## Installation

The pipeline can be set up in two easy steps:

> 1. Clone the github directory onto your system.

```
git clone https://github.com/Snitkin-Lab-Umich/snpkit.git

```

> 2. Use snpkit/environment_gubbins.yml files to create dependencies requried to run gubbins.

```
conda env create -f snpkit/envs/environment_gubbins.yml -n gubbins
```

Activate gubbins environment and install gubbins. Check gubbbins installation 

```
conda activate gubbins
conda install gubbins

run_gubbins.py --help
```

Activate snpkit conda environment. 
```
conda activate /nfs/turbo/umms-esnitkin/conda/snpkit

python snpkit/snpkit.py -h
```

## Quick Start

Lets say you want to detect variants for more than a few hundred samples against a reference genome KPNIH1 and want the pipeline to run in parallel on HPC cluster. 


- Run the first step of pipeline with option "-steps All" that will call variants for samples placed in test_readsdir against a reference genome KPNIH1

```

python snpkit/snpkit.py \
-type PE \
-readsdir /Path-To-Your/test_readsdir/ \
-outdir /Path/test_output_core/ \
-analysis output_prefix \
-index KPNIH1 \
-steps call \
-cluster cluster \
-scheduler SLURM \
-clean

```

- The above command will run the variant calling part of the pipeline on a set of PE reads residing in test_readsdir. 
- The results will be saved in the output directory test_output_core. 
- The reference genome and its path will be detected from the KPNIH1 settings that is set in config file.


The results of variant calling will be placed in an individual folder generated for each sample in the output directory. A log file for each sample will be generated and can be found in each sample folder inside the output directory. 

- Run the second part of the pipeline to generate SNP and Indel Matrices and various multiple sequence alignments outputs.

```
python snpkit/snpkit.py \
-type PE \
-readsdir /Path-To-Your/test_readsdir/ \
-outdir /Path/test_output_core/ \
-analysis output_prefix \
-index reference.fasta \
-steps parse \
-cluster cluster \
-scheduler SLURM \
-gubbins yes \ 
-mask 

```

This step will gather all the variant call results of the pipeline, generate SNP-Indel Matrices, qc reports and core/non-core sequence alignments that can be used as an input for phylogenetic analysis such as gubbins-iqtree-beast.

## Input

The pipeline requires three main inputs - readsdir, name of the reference genome and path to the config file.

**1. readsdir:** Place your Illumina SE/PE reads in a folder and give path to this folder with -readsdir argument. Apart from the standard Miseq/Hiseq fastq naming convention (R1_001_final.fastq.gz), other acceptable fastq extensions are: 

```

- R1.fastq.gz/_R1.fastq.gz, 
- 1_combine.fastq.gz, 
- 1_sequence.fastq.gz, 
- _forward.fastq.gz, 
- _1.fastq.gz/.1.fastq.gz.

```

**2. config:** A high level easy-to-write YAML format configuration file that lets you configure your system-wide runs and specify analysis parameters, the path to the installed tools, data and system-wide information.

- This config file will contain high-level information such as locations of installed programs like GATK, cores and memory usage for running on HPC compute cluster, path to a reference genome, various parameters used by different tools. These settings will apply across multiple runs and samples. 

- The config file stores data in KEY: VALUE pair. 

- An example [config](https://github.com/alipirani88/snpkit/blob/master/config) file with default parameters is included with the installation folder. You can customize this config file and provide it with the -config argument or edit this config file based on your requirements. 

- Parameters for each of the tools can be customized under the 'tool_parameter' attribute of each tool in config file. 

- If you wish to run the pipeline in a HPC compute environment such as PBS or SLURM, change the number of nodes/cores memory reuirements based on your needs else the pipeline will run with default settings.


**3. index:** a reference genome index name as specified in a config file. For example; if you have set the reference genome path in config file as shown below, then the required value for command line argument -index would be -index KPNIH1

```
[KPNIH1]
# path to the reference genome fasta file.
Ref_Path: /nfs/turbo/umms-esnitkin/data_sharing/reference/KPNIH1/
# Name of reference genome fasta file.
Ref_Name: KPNIH1.fasta
```

Here, Ref_Name is the reference genome fasta file located in Ref_Path. Similarly, if you want to use a different version of KPNIH1 reference genome, you can create a new section in your config file with a different index name.

```
[KPNIH1_V2024]
# path to the reference genome fasta file.
Ref_Path: /nfs/turbo/umms-esnitkin/data_sharing/reference/KPNIH1_V2024/
# Name of reference genome fasta file.
Ref_Name: KPNIH1_V2024.fasta
```

The pipeline assumes that you have placed the reference genome fasta file `KPNIH1.fasta` in folder `/nfs/turbo/umms-esnitkin/data_sharing/reference/KPNIH1/`, a genbank annotation file with extension `.gbf`, and PHASTEST phage regions results. The PHASTEST files that the pipeline expects are `summary.txt` and `phage_regions.fna`. Documentation on how to run PHASTEST on their web server or locally (Docker file) can be found here: https://phastest.ca/ 

Note: By default, Prokka outputs only .gbk files and not .gbf. You can change the extension of your genbank file to `.gbf` from `.gbk`. This could solve the possible extension requirement issue.

### For detailed information, please refer to the [wiki](https://github.com/Snitkin-Lab-Umich/snpkit/wiki) page.

