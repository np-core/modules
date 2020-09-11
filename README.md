# Nextflow Pipeline Modules (DSL2)

![](https://img.shields.io/badge/lang-dsl2-41ab5d.svg)
![](https://img.shields.io/badge/version-0.1.0-addd8e.svg)

Nextflow pipeline modules and workflows for nanopore-driven bacterial genomics in DSL2

Consider to configure your network setup for distributed pipeline operation and resourcing with [`Netflow`](https://github.com/np-core/netflow). Details for pipeline operations can be found in th corresponding repositories. 

## Pipelines

Purpose and usage of the main pipeline applications in `np-core`

### :peacock: np-signal

`nextflow run np-core/np-signal --help true`

Signal level analysis pipeline for Fast5 files, using NanoPath, Poremongo and Achilles linking into neural network training and adaptive sampling. Real-time basecalling and monitoring in the dashboard component Blobfish, the ugly deep sea cousin of MinKNOW.

### :drop_of_blood: np-pathogen

`nextflow run np-core/np-pathogen --help true`

Pipeline for detection and characteriziation of pathogens from metagenomic data, including real-time evaluation and reporting dashboard. Main component for the Queensland Genomics sepsis applications.

### :sauropod: np-phybeast

`nextflow run np-core/np-phybeast --help true`

Pipeline for phylogenomics and phylodynamics of bacterial pathogens using population-wide data, including subworkflows for outbreak attribution using nanopore data and Beastling implementation of bacterial Birth-Death model implementations on BEAGLE GPUs.

### :octopus: np-assembly

`nextflow run np-core/np-assembly --help true`

Pipeline for hybrid and nanopore  bacterial genome assembly including final genotyping and annotations of the assemblies. I use the pipeline for both Illumina population genomics (genotyping from large collection of genome assemblies) and to process hybrid sequence data from multiplex bacterial runs. 


## Modules

**Basecalling**

* Guppy
* qcat

**Quality control*

* nanoq
* rasusa
* checkm
* coverm
* fastp
* dnadiff

**Taxonomy**

* Kraken2
* Bracken
* Pathogen (links, regex, host removal, coverage)

**Variants**

* Medaka
* Megalodon
* Snippy
* Gubbins

**Assembly**

* Flye
* Raven
* Shovill
* Racon
* Pilon
* Medaka

**Phylodynamics**

* Phybeast (date randomisation, variant sites)
* RAxML-NG
* TreeTime

**Genotyping**

* Genotype (Mykrobe, Abricate, SCCion)
