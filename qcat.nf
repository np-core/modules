process Qcat {

    tag { id }
    label "ont"

    publishDir "$params.outdir/qcat", mode: "copy"

    input:
    tuple val(id), file(fq)

    output:
    tuple val(id), file("${id}")

    """
    qcat -f $fq -b ${id} --trim $params.qcat_params
    """

}

