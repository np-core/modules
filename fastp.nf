process Fastp {

    label "fastp"
    tag { id }

    input:
    tuple val(id), file(forward), file(reverse)

    output:
    tuple val(id), file("${id}_1_qc.fq.gz"), file("${id}_2_qc.fq.gz")

    """
    fastp -i $forward -I $reverse -o ${id}_1_qc.fq.gz -O ${id}_2_qc.fq.gz --thread $task.cpus
    """

}

process FastpTraining {

    label "fastp"
    tag { id }

    input:
    tuple val(model), val(id), file(forward), file(reverse), file(ont)

    output:
    tuple val(model), val(id), file("${id}_1_qc.fq.gz"), file("${id}_2_qc.fq.gz"), file(ont)

    """
    fastp --in1 $forward --in2 $reverse --out1 ${id}_1_qc.fq.gz --out2 ${id}_2_qc.fq.gz --thread $task.cpus
    """

}

process FastpEvaluation {

    label "fastp"
    tag { id }

    input:
    tuple val(eval_set), val(id), file(forward), file(reverse)

    output:
    tuple val(eval_set), val(id), file("${id}_1_qc.fq.gz"), file("${id}_2_qc.fq.gz")

    """
    fastp --in1 $forward --in2 $reverse --out1 ${id}_1_qc.fq.gz --out2 ${id}_2_qc.fq.gz --thread $task.cpus
    """

}