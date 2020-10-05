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

    publishDir "${params.outdir}/medaka", mode: "copy", pattern: "${id}.vcf"
    publishDir "${params.outdir}/medaka", mode: "copy", pattern: "${id}.txt"

    input:
    tuple val(id), file(bam), file(bai)
    file(reference)

    output:
    tuple val(id), file("${id}.vcf"), file("${id}.txt")
    tuple val(id), file(bam), file(bai)

    """
    medaka consensus --model $params.medaka_model --threads $task.cpus ${id}.bam ${id}.hdf
    medaka snp --threshold 1 $reference ${id}.hdf ${id}.vcf
    pysamstats -t variation_strand $bam -f $reference > ${id}.txt
    """

}

process MedakaVariantsTraining {

    label "medaka"
    tag { "$model_name" }

    publishDir "${params.outdir}/medaka/${model_name}/${reference.baseName}", mode: "copy", pattern: "*.vcf"
    publishDir "${params.outdir}/medaka/${model_name}/${reference.baseName}", mode: "copy", pattern: "*.txt"

    input:
    tuple val(model_name), val(coverage), val(ids), file(reference), file("bams/*"), file("bams/*")

    output:
    tuple val(model_name), file("*.vcf"), file("*.txt")

    """
    for i in $ids; do
        medaka consensus --model $params.medaka_model --threads $task.cpus bams/\${i}_${coverage}.bam \${i}_${coverage}.hdf
        medaka snp --threshold 1 $reference \${i}_${coverage}.hdf \${i}_${coverage}.vcf
        pysamstats -t variation_strand bams/\${i}_${coverage}.bam -f $reference > \${i}_${coverage}.txt
    done
    """

}
