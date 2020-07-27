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
* Phylodynamics

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
* Rraxml
* Shovill
* Snippy
* TreeTime
* Phybeast
