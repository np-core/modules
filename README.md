# Core Pipelines

![](https://img.shields.io/badge/lang-dsl2-41ab5d.svg)
![](https://img.shields.io/badge/version-0.1.0-addd8e.svg)

Nextflow pipeline modules and workflows for nanopore-driven bacterial genomics in DSL2

Before launching pipelines and dashboards please consider configuring your network setup for distributed resourcing with [`Netflow`](https://github.com/np-core/netflow)

## Pipelines

Purpose and usage of the main pipeline applications in `np-core`

### :octopus: Signal

`nextflow run np-core/phybeast --help`

Signal level analysis pipeline for Fast5 files, using NanoPath, Poremongo and Achilles linking into neural network training and adaptive sampling. Real-time basecalling and monitoring in the dashboard component Blobfish, the ugly deep sea cousin of MinKNOW.

### :dragon: Pathogen

Pipeline for detection and characteriziation of pathogens from metagenomic data, including real-time evaluation and reporting dashboard. Main component for the Queensland Genomics sepsis applications.

### :sauropod: Phybeast

`nextflow run np-core/phybeast --help`

Pipeline for phylogenomics and phylodynamics of bacterial pathogens using population-wide data, including subworkflows for outbreak attribution using nanopore data and Beastling implementation of bacterial Birth-Death model implementations on BEAGLE GPUs.

### :crocodile: Assembly

`nextflow run np-core/assembly --help`

Pipeline for hybrid and nanopore  bacterial genome assembly including final genotyping and annotations of the assemblies. I use the pipeline for both Illumina population genomics (genotyping from large collection of genome assemblies) and to process hybrid sequence data from multiplex bacterial runs. 


