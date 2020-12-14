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
    tuple val(id), file(fq)
    file(reference)

    output:
    tuple val(id), file("${id}.vcf"), file("${id}.txt")
    tuple val(id), file(bam), file(bai)

    """
    medaka_haploid_variant --model $params.medaka_model --threads $task.cpus --output_dir vars $fq $reference
    mv vars/consensus_to_ref.vcf ${id}_${coverage}.vcf
    pysamstats -t variation_strand vars/calls_to_draft.bam -f $reference > ${id}_${coverage}.txt
    """

}

process MedakaTraining {

    label "medaka"
    tag { "$model_name - $id - $reference" }

    memory { params.medaka_mem * task.attempt }

    errorStrategy { task.exitStatus in 1..143 ? 'retry' : 'ignore' }
    maxRetries 5

    publishDir "${params.outdir}/${ref}/polishers/variants", mode: "copy", pattern: "${id}_${coverage}.vcf"
    publishDir "${params.outdir}/${ref}/polishers/variants", mode: "copy", pattern: "${id}_${coverage}.txt"

    input:
    tuple val(model_name), val(id), val(ref), val(coverage), file(reference), file(fq), file(snippy_vcf)

    output:
    tuple val(model_name), val(ref), file("${id}_${coverage}.vcf"), file("${id}_${coverage}.txt"), file(snippy_vcf)


    """
    medaka_haploid_variant --model $params.medaka_model --threads $task.cpus --output_dir vars $fq $reference
    mv vars/consensus_to_ref.vcf ${id}_${coverage}.vcf
    pysamstats -t variation_strand vars/calls_to_draft.bam -f $reference > ${id}_${coverage}.txt
    """

}

process MedakaEvaluation {

    label "medaka"
    tag { "$id - $reference" }

    memory { params.medaka_mem * task.attempt }

    errorStrategy { task.exitStatus in 1..143 ? 'retry' : 'ignore' }
    maxRetries 5

    publishDir "${params.outdir}/${reference.simpleName}/evaluation/${eval_set}/medaka", mode: "copy", pattern: "${id}.vcf"
    publishDir "${params.outdir}/${reference.simpleName}/evaluation/${eval_set}/medaka", mode: "copy", pattern: "${id}.txt"

    input:
    tuple val(eval_set), val(id), file(fq)
    each file(reference)

    output:
    tuple val(eval_set), val(id), val("${reference.simpleName}"), file("${id}.vcf"), file("${id}.txt")


    """
    medaka_haploid_variant --model $params.medaka_model --threads $task.cpus --output_dir vars $fq $reference
    mv vars/consensus_to_ref.vcf ${id}.vcf
    pysamstats -t variation_strand vars/calls_to_draft.bam -f $reference > ${id}.txt
    """

}
