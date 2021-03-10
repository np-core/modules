process Prokka {

    tag { id }
    label "prokka"

    publishDir "$params.outdir/qcat", mode: "copy"

    input:
    tuple val(id), file(fasta)

    output:
    tuple val(id), file("${id}_prokka")

    """
    prokka --compliant --centre SLU --outdir ${id}_prokka \
        --locustag $locus_tag --prefix $locus_tag --kingdom Bacteria \
        --genus $params.genus --species $params.species --usegenus $fasta \

    """

}
