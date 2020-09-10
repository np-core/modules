# DSL2 pipelines

![](https://img.shields.io/badge/lang-dsl2-41ab5d.svg)
![](https://img.shields.io/badge/version-0.1.0-addd8e.svg)

Nextflow implementations as pipeline modules and workflows for nanopore-driven bacterial genomics in DSL2. Config files specify parameters and run profiles using `Singularity` containers described in the container repository and available to pull by corresponding name from [DockerHub](https://hub.docker.com/u/esteinig).


## Pipelines

Purpose and usage of the main pipeline applications in `np-core`

### ::octopus:: Signal

`nextflow run np-core/phybeast --help`

Signal level analysis pipeline for Fast5 files, using NanoPath, Poremongo and Achilles linking into neural network training and dashboard for adaptive sampling.

### :dragon: Pathogen

Pipeline for detection and characteriziation of pathogens fro memtagenomic data, inlcuding real-time evaluation and reporting dashboard in Vue. Main component for the Queenslan Genomics sepsis groups.

### :sauropod: Phybeast

`nextflow run np-core/phybeast --help`

Pipeline for phylogenomics and phylodynamics of bacterial pathogens using population-wide data, inlcuding subworkflow for outbreak attribution using nanopore data and Beastling implementation of bacterial Birth-Death model implementations on BEAGLE GPU.

### :crocodile: Assembly

`nextflow run np-core/assembly --help`

Pipeline for hybrid and nanopore  bacterial genome assembly including final genotyping and annotations of the assemblies. I use the pipeline for both Illumina population genomics (genotyping from large collection of genome assemblies) and to process hybrid sequence data from multiplex bacterial runs. 


