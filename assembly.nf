// Hybrid assembly workflow: DSL2

nextflow.enable.dsl=2

// Helper functions

def get_single_fastx( glob ){
    return channel.fromPath(glob) | map { file -> tuple(file.baseName, file) }
}

def get_paired_fastq( glob ){
    return channel.fromFilePairs(glob, flat: true)
}

def get_matching_data( channel1, channel2 ){

    // Get matching data by ID (first field) from two channels
    // by crossing, checking for matched ID and returning
    // a combined data tuple with a single ID (first field)

    channel1.cross(channel2).map { crossed ->
        if (crossed[0][0] == crossed[1][0]){
            tuple( crossed[0][0], *crossed[0][1..-1], *crossed[1][1..-1] )
        } else {
            null
        }
    }
    .filter { it != null }

}

// ONT Quality control subworkflow

params.fastq = 'fastq/*.fq'
params.outdir = 'PNG'
params.reference = "$PWD/jdk.fasta"
params.length = 200
params.quality = 7
params.subsample = 200
params.genome_size = '2.8m'
params.preset = 'minimap2-ont'  // coverm

ref = file(params.reference)

include { Nanoq } from './modules/nanoq'
include { Rasusa } from './modules/rasusa'
include { CoverM  } from './modules/coverm'
include { NanoqStatistics } from './modules/nanoq' addParams(length: 0, quality: 0)

workflow ont_qc {
    take:
        reads // id, reads
    main:
        NanoqStatistics(reads)
        Nanoq(reads)
        CoverM(Nanoq.out, ref)
        Rasusa (Nanoq.out) 
    emit:
        Rasusa.out // id, reads
}

// ONT assembly and polishing subworkflow

params.assembly_options = "--plasmids"
params.medaka_model = "r941_min_high_g360"

include { Flye } from './modules/flye'
include { Racon } from './modules/racon'
include { Medaka } from './modules/medaka'

workflow ont_assembly {
    take:
        reads // id, reads
    main:
        Flye(reads)
        get_matching_data(Flye.out, reads) | Racon
        Medaka(Racon.out)
    emit:
        Flye.out  // id, assembly
        Medaka.out  // id, polished assembly
}

// Illumina assembly, hybrid correction and reference comparison subworkflow

params.fasta = "fasta/*.fasta"
params.illumina = "fastq/*_{1,2}.fastq.gz"
params.depth = 200
params.assembler = "skesa"

params.tag = null
params.saureus = true
params.kpneumoniae = false

include { Fastp } from './modules/fastp'
include { Shovill } from './modules/shovill'

// Assign tags for separate output folders in params.outdir / hybrid / params.tag
include { Pilon as FlyePilon } from './modules/pilon' addParams( tag: 'flye' )
include { Pilon as MedakaPilon } from './modules/pilon' addParams( tag: 'medaka' )

include { Dnadiff as MedakaComparison } from './modules/dnadiff' addParams( tag: 'medaka' )
include { Dnadiff as FlyeHybridComparison } from './modules/dnadiff' addParams( tag: 'flye_hybrid' )
include { Dnadiff as MedakaHybridComparison } from './modules/dnadiff' addParams( tag: 'medaka_hybrid' )

include { Genotype as IlluminaGenotype } from './modules/genotype' addParams( tag: 'illumina' )
include { Genotype as AssemblyGenotype } from './modules/genotype' addParams( tag: 'preassembled' )
include { Genotype as MedakaHybridGenotype } from './modules/genotype' addParams( tag: 'hybrid' )
include { Genotype as MedakaGenotype } from './modules/genotype' addParams( tag: 'ont' )


workflow hybrid_genotype_assembly {
   // Polished ONT assembly and Illumina reference assembly
   get_single_fastx(params.fastq) | ont_qc | ont_assembly
   get_single_fastx(params.fasta)  | AssemblyGenotype // from preassembled genomes
   
   println ont_assembly.out[1] | MedakaGenotype  // access second emission medaka of the ont-qc workflow

   get_paired_fastq(params.illumina) | Fastp | Shovill | IlluminaGenotype  // illumina reference assembly and genotypes

   // Branch into correction and comparison to Illumina from Flye and Medaka
   get_matching_data(ont_assembly.out[0], Fastp.out) | FlyePilon
   get_matching_data(ont_assembly.out[1], Fastp.out) | MedakaPilon | MedakaHybridGenotype
   
   // Compare Flye and Medaka assemblies to Illumina reference assembly
   get_matching_data(ont_assembly.out[1], Shovill.out) | MedakaComparison
   get_matching_data(FlyePilon.out, Shovill.out) | FlyeHybridComparison
   get_matching_data(MedakaPilon.out, Shovill.out) | MedakaHybridComparison
}

// Execute

workflow {
    hybrid_genotype_assembly()
}