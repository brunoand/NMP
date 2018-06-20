# NMP

NMP is a pipeline for Microbiome Quality check, trimming and Taxonomic annotation.

Usage = nextflow run -c nextflow.config NMP.nf --reads_R1 R1 --reads_R2 R2 --date date --prefix mysample --outdir path -with-docker Docker_image

nextflow required version 0.29.0.4803

NMP uses the software bbDUK for trimming and decontamination (host genome) steps.
The latest build of Metaphlan2 is used for the taxonomic classification step
FastQC is used to retrive quality metrics of fastqs (Before and after the quality check step)

All results are stored in 3 folders inside the outdir/prefix

Matrix - where all the matrices of a given pattient is stored as $date.tsv (if they have the same preffix)
QC - where all the Fastqc files and trimmed fastq files are stored
Plots - where barplots and their legends are stored

All softwares needed for the analysis are stored in a Docker image

This pipeline relies on two different files, NMP.nf and nextflow.config:

NMP.nf is the pipeline per se, all scripts and processes are stored inside this file.
nextflow.config is the configuration file. This file gives you control over all the softwares used by this pipeline, such as:

External files pathways (adapters, artifacts, Reference genome (indexed)
Software parameters, like minimum quality for read trimming, bowtie2 options and others

And more important, you can define the maximum time, memory and number of cores for each proccess of the pipeline.
