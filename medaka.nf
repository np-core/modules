process Medaka {

    tag { id }
    label "medaka"

    publishDir "$params.outdir/ont/assembly", mode: "copy", pattern: "*.medaka.fasta"

    input:
    tuple val(id), file(racon_assembly), file(fastq)

    output:
    tuple val(id), file("${id}.medaka.fasta")
    
    """ 
    medaka consensus -i $fastq -d $racon_assembly -o racon_medaka -t $task.cpus -m $params.medaka_model
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
    tag { "$model_name - $id - $reference" }

    publishDir "${params.outdir}/medaka/${model_name}/${reference.baseName}", mode: "copy", pattern: "${id}_${coverage}.vcf"
    publishDir "${params.outdir}/medaka/${model_name}/${reference.baseName}", mode: "copy", pattern: "${id}_${coverage}.txt"

    input:
    tuple val(model_name), val(id), val(coverage), file(reference), file(bam), file(bai), file(snippy_vcf)

    output:
    tuple val(model_name), val("${reference.baseName}"), file("${id}_${coverage}.vcf"), file("${id}_${coverage}.txt"), file(snippy_vcf)

    """
    medaka consensus --model $params.medaka_model --threads $task.cpus $bam ${id}_${coverage}.hdf
    medaka snp --threshold 1 $reference ${id}_${coverage}.hdf ${id}_${coverage}.vcf
    pysamstats -t variation_strand $bam -f $reference > ${id}_${coverage}.txt
    """

}
