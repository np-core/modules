process GraftM {

    label "graftm"
    tag { "$id" }

    publishDir "${params.outdir}/graftm/$pkg", mode: "copy", pattern: "*"

    input:
    tuple val(id), file(fwd), file(rev), val(pkg), file(package)

    output:
    file("${id}")

    """
    graftm graft --forward $fwd --reverse $rev --graftm_package $package --output_directory ${id}
    """

}