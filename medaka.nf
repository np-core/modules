process Medaka {

    tag { id }
    label "medaka"

    publishDir "$params.outdir/ont/assembly", mode: "copy", pattern: "*.medaka.fasta"

    input:
    tuple val(id), file(racon_assembly), file(fastq)

    output:
    tuple val(id), file("${id}.medaka.fasta")
    
    """ 
    medaka_consensus -i $fastq -d $racon_assembly -o racon_medaka -t $task.cpus -m $params.medaka_model
    mv racon_medaka/consensus.fasta ${id}.medaka.fasta
    """

}

process MedakaVariants {

    label "medaka"
    tag { "$id" }

    publishDir "${params.outdir}/medaka", mode: "copy"

    input:
    tuple val(id), file(fastq)

    output:
    tuple val(id), file("${id}.vcf")
    file("${id}.bam")

    """
    minimap2 -ax map-ont $reference $fastq | samtools sort | samtools view -Sb > ${id}.bam
    samtools index ${id}.bam
    medaka consensus --model $params.medaka_model --threads $task.cpus ${id}.bam ${id}.hdf
    medaka snp --threshold 1 $reference ${id}.hdf ${id}.vcf
    """

}
