process CheckM {

    tag { id }
    label "coverm"

    publishDir "$params.outdir/ont/qc", mode: "copy", pattern: "*.txt"

    input:
    tuple val(id), file(fa)
    
    output:
    tuple val(id), file("${id}.checkm.txt")

    """
    checkm lineage_wf -x fa $PWD checkm_output
    mv checkm_output/
    """

}