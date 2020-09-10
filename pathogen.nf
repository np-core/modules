#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
==============================================================================================================================
                                        N P - P A T H O G E N   P I P E L I N E
==============================================================================================================================

 Metagenomic pathogen detection and quality control for online or offline
 ONT and Illumina (DNA or RNA) sequencing data

 Documentation: https://github.com/np-core/pathogen

 Original development by Queensland Genomics, Australian Institute of Tropical Health 
 and Medicince, The Peter Doherty Intitute for Infection and Immunity

Developers:

    Eike Steinig  @esteinig  < @EikeSteinig >  ont, dna, online, bacterial, clinical, reports, netflow, nanopath
    Cadhla Firth  @cfirth    < @cfirth      >  ont, rna, illumina, viral, metatransciptomics, clinical, lab workflows
    Lachlan Coin  @ljcoin    < @ljmcoin     >  ont, dna, rna, illumina, hybrid, metagenomics, online, clinical, reports

Contributors:   

    Tania Duarte        UQ      @tduarte        < @TaniaDuarte >  workflow, sputum labwork, cystic fibrosis applications
    Ross Balch          UQ      @rbalch         < @RossBalch   >  workflow, blood labwork, sepsis applications
    Masalgi Mills       JCU     @mmills         < @mmills      >  database, human global genome project assembly
    Matt Field          JCU     @mattfield      < @FieldM      >  workflow, coverage mapping component quality control
    Daniel Rawlinson    PDIII   @dna-raw        < @DRawlinson  >  worfklow, container management 

Pipeline part of the NanoPath core framework:

    https://github.com/np-core

NanoPath distributed pipeline framework Netflow:

    https://github.com/np-core/netflow

For interactive dashboard operation and report generation:

    https://github.com/np-core/nanopath

----------------------------------------------------------------------------------------
*/

import java.nio.file.Paths;

nextflow.enable.dsl=2;

// nextflow config is loaded before compilation

params.files = "$PWD"  // input file directory containing ONT (.fq, .fastq) or Illumina (.fq.gz, .fastq.gz, PE) reads
params.recursive = false // recursive file search
params.outdir = "$PWD/results" // result output directory, must be full path!

params.online = false  // online (approximate and continous Raven assembly) or offline (complete metaFlye assembly) operations
params.read_regex = /barcode\d\d/ // group individual input files by regex (online) and execute pipelines for each group (e.g. barcode)

params.min_length = 50  // minimum read length filter to remove smashed uninformative fragments (mild)
params.min_quality = 5  // minimum quality length filter to remove low quality fragments (mild)

params.assembly_every = 2 // assemble the complete read collection for each group every <assembly_every> batches incoming from online sequencing run (usually 4000 reads) - host removed in batch before
params.assembly_update = "${params.outdir}/fastq"  // this is a bit unorthodox, but symlinks fastqs into this dir, to compute repeated mag assembly on
params.assembly_options = "" // metagenome assembler options for either raven (ont dna online) metaflye (ont dna offline)
params.assembly_markers = ['vfdb', 'resfinder']  // abricate databases to search for marker genes on assembled contigs / MAGs

params.host_database = "$HOME/resources/host" // kraken database to use for deep host cleaning, only read classificaion not abundance
params.host_reference = "$HOME/resources/host.ont.mmi"  // minimap2 index to use for deep host cleaning, via additional read mapping

params.taxonomy_databases = ["$HOME/resources/minikraken2"] // taxonomy databases to use for classification
params.coverage_databases = ["$HOME/resources/refseq.fasta"]  // database of reference genomes to use for coverage computation

params.bracken_threshold = 5  // minimum number of reads to be assigned to level to receive reads in abundance computation
params.bracken_level = "S"  // level of taxonomic classification that wil lreceive reads from other levels for abundance estimation
params.bracken_length = 100  // read length for kmer indices used in bracken, files must be present with the databases

params.coverage_threshold = 1000  // threshold of bracken abundance reads to initiate coverage mapping, prevents coverage estimates for everything detected in complex samples, ref genome selection automated
params.coverage_window = 500  // length of coverage window to bin coverage estimates along reference genome
params.coverage_options = "--mapq 60" // mosdepth settings, e.g. to gilter minimap2 alignments before coverage estimate

params.rna = false
params.illumina = false
params.tail = "_{1,2}"


// Workflow version

version = '0.1.7'

def helpMessage() {

    log.info"""
    =========================================
    Metagenomic pathogen detection v${version}
    =========================================

    Usage:

    The typical command for running the pipeline is as follows for real-time (batch-wise) nanopore processing (DNA):

        nextflow run np-core/pathogen --files fastq_input_directory --host_db /db/hominid --databases /db/minikaken2

    Input / output:

        --files                 the path to directory of fastq files: illumina pe (.fastq, .fastq.gz) ont (.fq, .fq.fz) [${params.files}]
        --recursive             activate recrusive file search for input directory [${params.recursive}]
        --outdir                the path to directory for result outputs [${params.outdir}]

    Nanopore DNA read file groupings (produced from basecaller in batches, usually labeled with barcode) and read filters:

        --online                ont dna, ont rna: activate online processing for batch-wise read files (default 4000 reads from MinKNOW) [${params.online}]
        --read_regex            ont dna, ont rna: group files by filename regex e.g. group by barcode [${params.read_regex}]

        --min_length            ont dna, ont rna: filter reads by minimum length before processing [${params.min_length}]
        --min_quality           ont dna, ont rna: filter reads by minimum quality before processing [${params.min_quality}]

    Database arguments (comma delineated strings for lists on command-line):

        --host_database         the specific host database to use for deep cleaning reads (Kraken2) [${params.host_database}]
        --host_reference        the specific host reference genome index to use for deep cleaning reads (minimap2) [${params.host_reference}]

        --taxonomy_databases    list of taxonomic databases to employ for read or transcript classification (Kraken2) ${params.taxonomy_databases}
        --coverage_databases    list of reference genome databases to query when assessing genome coverage (NanoPath) ${params.coverage_databases}

    Taxonomy abundance arguments:

        --bracken_level         ont dna, illumina dna: compute abundance adjusted read counts at this taxonomic level [${params.bracken_level}]
        --bracken_length        ont dna, illumina dna: bracken read length parameter, kmer indices must be present with databases [${params.bracken_length}]
        --bracken_threshold     ont dna, illumina dna: minimum number of reads required to be distributed reads in abundance estimate [${params.bracken_threshold}]

    Metagenome assembly arguments:

        --assembly_options      ont dna (online, offline): arguments to pass to workflow metagenome assembler [${params.assembly_options}]
        --assembly_markers      ont dna (online, offline): scan the metagenome assembled contigs for marger genes [${params.assembly_markers}]

        --assembly_update       ont dna (online): assemble the entire metagenome from fastq files available in this path (symlinked in pipeline) [${params.assembly_update}]
        --assembly_every        ont dna (online): assemble the entire metagenome every <n> batches during real-time processing, using fast assembly with Raven [${params.assembly_every}]

    Coverage mapping arguments:
    
        --coverage_threshold    ont dna, illumina dna: bracken abundance threshold in reads by bracken to trigger reference coverage mapping [${params.coverage_threshold}]
        --coverage_window       ont dna, illumina dna, illumina rna: estimate genome coverage across the reference genome windows of this size [${params.coverage_window}]
        --coverage_options      ont dna, illumina dna, illumina rna: setting to mosdepth for example to filter reads by mapping quality ["${params.coverage_options}""]

    Pipeline subworkflow arguments:

        --rna                   activate rna sequencing pipeline for illumina or ont pipelines [${params.rna}]
        --illumina              search for illumina files (.fastq.gz, .fq.gz) in path and activate short-read pipeline [${params.illumina}]
        --tail                  glob the illumina fastq input files by their paired end tails ["${params.tail}""]

    """.stripIndent()
}


params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Database file checks, if not exist exit with error:

def check_file(file) {
    
    path = Paths.get(file)

    if (path.exists()){
        log.info"""
        Detected input file: $file
        """
    } else {
        log.info"""
        Failed to detect input file: $file
        """
        exit 0
    }
}


// Important for real-time assembly!
// Create the output directory and the assembly update directory:

File assembly_update_directory = new File(params.assembly_update)
assembly_update_directory.mkdirs() // Create the update dir so it ...
assembly_update = file(params.assembly_update) // ... can be declared as input


check_file(params.host_database)
check_file(params.host_reference)

params.taxonomy_databases.each { check_file(it) }
params.coverage_databases.each { check_file(it) }

// Parameter checks and variable configuratons

host_database = file(params.host_database)
host_reference = file(params.host_reference)

if ( params.taxonomy_databases instanceof String ){
    db = params.taxonomy_databases.split(",").collect { file(it) }
} else {
    db = params.taxonomy_databases.collect { file(it) }
}

if ( params.coverage_databases instanceof String ){
    refdb = params.coverage_databases.split(",").collect { file(it) }
} else {
    refdb = params.coverage_databases.collect { file(it) }
}


// Input file checks, if none detected exit with error:


if (params.recursive){
    _fastq_ont = ["${params.files}/**/*.fq", "${params.files}/**/*.fastq"]
    _fastq_illumina = ["${params.files}/**/*${params.tail}.fq.gz", "${params.files}/**/*${params.tail}.fastq.gz"]
} else {
    _fastq_ont = ["${params.files}/*.fq", "${params.files}/*.fastq"]
    _fastq_illumina = ["${params.files}/*${params.tail}.fq.gz", "${params.files}/*${params.tail}.fastq.gz"]
}



// Helper functions

def get_batch_fastq(glob){
    // For Raven assembly by read batches
    batch = 0
    return channel.fromPath(glob) | map { batch += 1; tuple(it.baseName, it, batch) }
}
def get_single_fastq(glob){
    return channel.fromPath(glob) | map { file -> tuple(file.baseName, file) }
}
def get_paired_fastq(glob){
    return channel.fromFilePairs(glob, flat: true)
}

include { LinkAssemblyOnline } from './modules/pathogen'
include { NanoqOnline } from './modules/nanoq'
include { DeepHostRemovalOnline } from './modules/pathogen'
include { RavenRegexOnline } from './modules/raven' 
include { KrakenOnline as ReadKrakenOnline } from './modules/kraken'
include { KrakenAssemblyOnline as RavenKrakenOnline } from './modules/kraken'

// Online MAG assembly evey few batches

workflow ont_assembly_online {
    take:
        read_batch
    main:   
        // Continous metagenome assembly across barcode regex
        LinkAssemblyOnline(read_batch)  // dummy process to symlink fastq into output directory for updated mag assembly < params.assembly_update > - ensures presence of files before metagenome assembly
        RavenRegexOnline(LinkAssemblyOnline.out, assembly_update) // fast mag assembly with raven by barcode every < params.assembly_every > batches, across all regex groups (barcoded) in linked fastq
        RavenKrakenOnline(RavenRegexOnline.out, db)  // contig classification for approximate mags
    emit:
        RavenKrakenOnline.out
}


// Online (batch-wise) classification 

workflow ont_dna_online {
    take:
        read_batch // id, fq, batch #
    main:
        NanoqOnline(read_batch)
        DeepHostRemovalOnline(NanoqOnline.out[0], host_database, host_reference) // hard filter on read quality > 5 and read length > 50 bp to remove highly fragmented reads (if called at all)
        ReadKrakenOnline(DeepHostRemovalOnline.out[0], db) // taxonomic classification of reads and abundance estimation, kraken + bracken 
    emit:
        ReadKrakenOnline.out
}


// Workflow to launch assembly and coverage mapping from dashboard

workflow ont_dna_launch {
    take:
        reads // id, fq
    main: 
        ReadKrakenOffline(reads, db) // taxonomic classification of reads and abundance estimation, kraken + bracken
        // fast or high quality metagenome assembly
        if (params.assembler == "raven"){
            contigs = MetaRavenRegex()
        } else {
            contigs = MetaFlyeRegex(DeepHostRemovalOffline.out[0])
        }
        FlyeKrakenOffline(contigs[0], db)  // mag taxonomic classification 
    emit:
        MetaFlyeOffline.out[0]
        MetaFlyeOffline.out[1]
        ReadKrakenOffline.out
        FlyeKrakenOffline.out
}



// Offline (non-batch, completed) ONT DNA

include { DeepHostRemovalOffline } from './modules/pathogen'
include { MetaFlye as MetaFlyeOffline } from './modules/flye'
include { NanoqOffline } from './modules/nanoq'
include { Kraken as ReadKrakenOffline } from './modules/kraken'
include { Kraken as FlyeKrakenOffline } from './modules/kraken'
include { FileGroupRegex } from './modules/pathogen'

workflow ont_dna_offline {
    take:
        reads // id, fq
    main: 
        NanoqOffline(reads) // emitted / analysed individually
        DeepHostRemovalOffline(NanoqOffline.out[0], host_database, host_reference) // remove host reads deep clean with kraken host and reference mapping
        ReadKrakenOffline(DeepHostRemovalOffline.out[0], db) // taxonomic classification of reads and abundance estimation, kraken + bracken
        MetaFlyeOffline(DeepHostRemovalOffline.out[0])  // high quality metagenome assembly
        FlyeKrakenOffline(MetaFlyeOffline.out[0], db)  // mag taxonomic classification 
    emit:
        MetaFlyeOffline.out[0]
        MetaFlyeOffline.out[1]
        ReadKrakenOffline.out
        FlyeKrakenOffline.out
}

workflow {
    if (params.illumina){
        if (params.rna){
            get_paired_fastq(_fastq_illumina) | illumina_dna // transcript abundance and taxonomic transcript classification, coverage mapping
        } else {
            get_paired_fastq(_fastq_illumina) | illumina_rna  // taxonomic abundance and mag assembly, coverage mapping
        }
    } else {
        if (params.rna){
            get_single_files(_fastq_ont)  // not implemented yet: ONT RNA
        } else {
            if (params.online){
                get_batch_fastq(_fastq_ont) |  ont_dna_online                                 // online barcoded taxonomic abundance and real-time (approximate) mag assembly, coverage mapping
            } else {
                get_single_fastq(_fastq_ont) |  ont_dna_offline                               // offline barcoded taxonomic abundance and approximate mag assembly, coverage mapping
            }
        }
    }
}



// // MAG binning approach

// include { Raven as BinAssembly } from './modules/raven' 
// include { KrakenONT as BinKraken  } from './modules/kraken'
// include { CheckM as BinCheckM } from './modules/checkm'

// workflow fast_mag_binning {
//     take:
//         assembly
//     main:
//         bins = Metabat2(assembly.out) | flatMap  // binning contigs across metagenome
//         BinAssembly(bins) // assemble individual bins    
//         BinKraken(BinRavenAssembly.out) // tax assignment  bin and assembly quality
//         BinCheckM(BinRavenAssembly.out) // mag quality of bins and assemblies
//     emit:
//         BinKraken.out
//         BinCheckM.out
// }

// include { SelectReference as SelectReferenceONT } from './modules/pathogen'
// include { SelectReference as SelectReferenceIllumina } from './modules/pathogen'
// include { SingleFastqLinkage } from './modules/pathogen'

// include { Fastp } from './modules/fastp'
// include { KrakenONT } from './modules/kraken'
// include { KrakenIllumina as KrakenIlluminaReads } from './modules/kraken'
// include { KrakenIllumina as KrakenIlluminaBins } from './modules/kraken'

// include { CoverageWindowONT } from './modules/pathogen'
// include { CoverageWindowIllumina as CoverageWindowIlluminaReadReference} from './modules/pathogen'
// include { CoverageWindowIllumina as CoverageWindowIlluminaBins } from './modules/pathogen'

// workflow kraken_ont_dna {
//     take:
//         reads // id, fq
//     main:
//         SingleFastqLinkage(reads)
//         RemoveHost(SingleFastqLinkage)
//         KrakenONT(RemoveHost.out, db) // Kraken2 + Bracken
//         ref = SelectReferenceONT(KrakenONT.out) | flatMap  // select references for taxa > params.coverage_threshold (reads)
//         CoverageWindowONT(reads, ref)  // online coverage mapping against individual selected references
//     emit:
//         KrakenONT.out, // id, kout, kreport, bout, breport
//         CoverageWindowONT.out // continously updated coverage file for all taxa over read threshold
// }     

// // Illumina PE DNA and RNA

// // Metagenome assembly pipeline similar to metaWRAP (but more control over parameters)

// workflow meta_assembly_illumina {
//     take:
//         reads
//     main:
//         MetaAssembly(reads)
//         Binning(MetaAssembly.out)
//         BinRefinement(Binning.out)
//         BinAssembly(BinRefinement.out)
//         BinAbundance(BinAssembly.out)
//     emit:
//         BinRefinement.out,
//         BinAssembly.out,
//         BinAbundance.out
// }

// workflow kraken_illumina_dna {
//     take:
//         reads // id, fwd, rev
//     main:
//         Fastp(reads)
//         if (params.illumina_kraken){
//             // Analogous to mapping approach in real-time nanopore for taxon coverage
//             KrakenIlluminaReads(Fastp.out, db) // output to params.outdir/kraken/{db}
//             ref = SelectReferenceIllumina(KrakenIllumina.out) | flatMap 
//             CoverageWindowIlluminaReads(reads, ref) // output to params.outdir/kraken/{db}/coverage
//         }
//         if (params.illumina_assembly){
//            meta_assembly_illumina(Fastp.out)  //output to
//         }  
//     emit:
        
// }    


// Abundance problem because expression varies in metatransciptome
// Basically does metagenome binning and then reassembly with metaWRAP
// submodules. 

// workflow kraken_illumina_rna {
//     take:
//         reads // id, fwd, rev
//     main:
        
//     emit:

// }    