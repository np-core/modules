process Dnadiff {

    tag { id }
    label "dnadiff"

    publishDir "$params.outdir/${params.tag}/dnadiff", mode: "copy", pattern: "*.report"

    input:
    tuple val(id), file(query_assembly), file(reference_assembly)

    output:
    tuple val(id), file("${id}.${params.tag}.report")


    """
    dnadiff $reference_assembly $query_assembly -p hybrid
    mv hybrid.report ${id}.${params.tag}.report
    """

}