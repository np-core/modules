# DSL2 pipelines

![](https://img.shields.io/badge/lang-dsl2-41ab5d.svg)
![](https://img.shields.io/badge/version-0.1.0-addd8e.svg)

Nextflow implementations as pipeline modules and workflows for nanopore-driven bacterial genomics in DSL2. Config files specify parameters and run profiles using `Singularity` containers described in the container repository and available to pull by corresponding name from [DockerHub](https://hub.docker.com/u/esteinig).

Pipelines:

* [assembly](https://github.com/np-core/assembly)
* [phybeast](https://github.com/np-core/phybeast)
* [metagenomics](https://github.com/np-core/metagenomics)
* [signal](https://github.com/np-core/signal)
* [sepsis](https://github.com/np-core/sepsis)
* [sketchy](https://github.com/np-core/sketchy)

Subworkflows:

* ONT QC
* ONT Assembly
* Illumina Assembly
* Reference Assembly
* Lineage Phylodynamics
* ONT Outbreak Attribution

Modules:

* CoverM
* dnadiff
* Fastp
* Flye
* Genotype
* Gubbins
* Medaka
* nanoq
* Pilon
* Racon
* Rasusa
* RAxML-NG
* Shovill
* Snippy
* TreeTime
* Phybeast

## Pipelines

Purpose and usage of the main pipeline applications in `np-core`:

### :crocodile: Assembly

`nextflow run np-core/assembly --help`

Pipeline for hybrid and nanopore  bacterial genome assembly including final genotyping and annotations of the assemblies. Usually run on barcoded sequence runs which have been basecalled and demultiplexed. I use the pipeline for both Illumina population genomics (genotyping from large collection of genome assemblies) and to process hybrid sequence data from multiplex bacterial runs. Illumina data is used to polish the long-read assemblies with Medaka.

### :sauropod: Phybeast

`nextflow run np-core/phybeast --help`

Pipeline for phylogenomics and phylodynamics of bacterial pathogens using population-wide data, inlcuding subworkflow for outbreak attribution using nanopore data alone.
