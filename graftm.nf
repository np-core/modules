process GraftM {

    label "graftm"
    tag { "$id" }

    publishDir "${params.outdir}/graftm", mode: "symlink", pattern: ""

    input:
    tuple val(id), file(fwd), file(rev)

    """
    graftm
    """

}