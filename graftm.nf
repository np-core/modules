process GraftM {

    label "graftm"
    tag { "$id" }

    publishDir "${params.outdir}/graftm/$pkg", mode: "copy", pattern: "*"

    input:
    tuple val(id), file(fwd), file(rev), val(pkg), file(graftm)

    output:
    file("${id}")

    """
    graftM graft --forward $fwd --reverse $rev --graftm_package $graftm --output_directory $id
    """

}

process GraftMAG {

    label "graftm"
    tag { "$id" }

    publishDir "${params.outdir}/graftm/$pkg", mode: "copy", pattern: "*"

    input:
    tuple val(id), file(fa), val(pkg), file(graftm)

    output:
    file("${id}")

    """
    graftM graft --forward $fa --graftm_package $graftm --output_directory $id
    """

}