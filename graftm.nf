process GraftM {

    label "graftm"
    tag { "$id" }

    publishDir "${params.outdir}/graftm", mode: "copy", pattern: "*"

    input:
    tuple val(id), file(fwd), file(rev)
    each file(package)

    output:
    file("${id}_$package_name")

    script:

    package_name = package.getName()

    """
    graftm graft --forward $fwd --reverse $rev --graftm_package $package --output_directory ${id}_$package_name
    """

}