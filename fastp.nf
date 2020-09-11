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