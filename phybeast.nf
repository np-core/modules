// Phybeast: bacterial phylodynamics workflow: DSL2

nextflow.enable.dsl=2

// Helper functions

def get_single_file( glob ){
    return channel.fromPath(glob) | map { file -> tuple(file.baseName, file) }
}

def get_paired_files( glob ){
    return channel.fromFilePairs(glob, flat: true)
}

// Non-recombinant core-genome SNP alignment


// Declare external files
params.reference = "$PWD/jdk.fasta"
reference = file(params.reference)

params.meta_data = "$PWD/meta_data.tsv"
meta_data = file(params.meta_data)

params.outdir = "$PWD/test_out"
params.tree_model = "GTR+G+ASC_LEWIS"

params.replicates = 100 // date randomisation test

include { Fastp } from './modules/fastp'
include { SnippyFastq } from './modules/snippy'
include { SnippyFasta } from './modules/snippy'
include { SnippyCore  } from './modules/snippy'
include { Gubbins  } from './modules/gubbins'

workflow snippy_fastq {
    take:
        reads // id, forward, reverse
    main:
        Fastp(reads)
        SnippyFastq(Fastp.out, reference)
    emit:
        SnippyFastq.out // id, results
}       

workflow snippy_fasta {
    take:
        contigs // id, fasta
    main:
        SnippyFasta(contigs, reference)
    emit:
        SnippyFasta.out // id, results
}  

workflow snippy_core {
    take:
        snippy // results
    main:
        SnippyCore(snippy.collect(), reference)
        Gubbins(SnippyCore.out) // wgs snp alignment
    emit:
        Gubbins.out // non-recombinant snp core alignment
}

include { RAxML } from './modules/raxml'
include { TreeTime } from './modules/treetime'
include { VariantSites } from './modules/phybeast'
include { DateRandomisation } from './modules/phybeast'

// Basic phylogeny and phylodynamics based on ML

workflow phylodynamics_ml {
    take:
        alignment
    main:
        VariantSites(alignment)
        RAxML(VariantSites.out)
        TreeTime(RAxML.out, meta_data, alignment)
        DateRandomisation(RAxML.out, TreeTime.out[0], meta_data, alignment)
    emit:
        RAxML.out
}

// Advanced models on GPU using BEAST2 and BEAGLE

// include { BirthDeathSkyline } from './modules/beastling'
// include { CoalescentSkyline } from './modules/beastling'
// include { MultiTypeBirthDeath } from './modules/beastling'

// // Should be used for exploratory runs unless sure that priors are configured appropriately

// params.iterations = 1000000  
// params.coupled_mcmc = false
// params.beast_options = "--beagle_gpu"

// params.cosky_config = "default"
// params.cosky_dimensions = [2, 4, 8, 16]

// params.bdss_config = "default"
// params.bdss_dimensions = [5, 6, 7, 8]

// params.mtdb_config = "default"
// params.mtbd_types = ['mssa', 'clade']

// include { Beast as BeastCosky } from './modules/beast'
// include { Beast as BeastBDSS } from './modules/beast'
// include { Beast as BeastMTBD } from './modules/beast'

// workflow phylodynamics_beast {
//     take:
//         core_snp_alignment
//     main:
//         CoalescentSkyline() | BeastCosky
//         BirthDeathSkyline() | BeastBDSS
//         MultiTypeBirthDeath() | BeastMTBD
// }

// Execute

params.subworkflow = "megalodon"
params.panels = "$HOME/LINEAGES/ST93/Megalodon"
params.candidates = "$HOME/LINEAGES/ST93/core.vcf"
params.devices = "1"
params.guppy_server_path = "/opt-guppy/bin/guppy_basecall_server"
params.guppy_params = "--trim_barcodes --chunk_size 512 --chunks_per_runner 2048 --gpu_runners_per_device 4"
params.guppy_config = "dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg"

include { MegalodonVariants } from './modules/megalodon'

def get_barcode_fast5(dir){

    port = 5555
    return channel.fromPath("$params.panels/**/*", type: 'dir').map { port += 1; tuple(port, it.getParent().getName(), it.getName(), it) }
}

        
workflow megalodon_variants {
    take:
        barcodes
    main:
        MegalodonVariants(reference, barcodes)
    emit:
        MegalodonVariants.out
}

workflow {
    if (params.subworkflow == "assembly"){

        fasta = get_single_file("FNQ*.fasta") | snippy_fasta
        fastq = get_paired_files("*_{R1,R2}.fastq.gz") | snippy_fastq
        fasta.mix(fastq) | snippy_core | phylodynamics_ml
    
    } else if (params.subworkflow == "megalodon"){

        get_barcode_fast5(params.panels) | megalodon_variants

    }
}