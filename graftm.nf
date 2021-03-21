process GraftM {

    label "graftm"
    tag { "$id" }

    publishDir "${params.outdir}/graftm", mode: "copy", pattern: "*"

    input:
    tuple val(id), file(fwd), file(rev)
    file(package)

    """
    graftm graft --forward $fwd --reverse $rev --graftm_package $package --output_directory ${id}
    """

}