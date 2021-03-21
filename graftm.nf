process GraftM {

    label "graftm"
    tag { "$id" }

    publishDir "${params.outdir}/graftm/$pkg", mode: "copy", pattern: "*"

    input:
    tuple val(id), file(fwd), file(rev), file(graftm), val(pkg)

    output:
    file("${id}")

    """
    graftm graft --forward $fwd --reverse $rev --graftm_package $graftm --output_directory $id
    """

}