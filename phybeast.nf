process DateRandomisation {

    label "phybeast"
    tag { "$tree" }

    publishDir "${params.outdir}/date_randomisation", mode: "copy"

    input:
    tuple val(id), file(tree)
    file(rate)
    file(meta_data)
    file(alignment)

    output:
    file("${id}_rates.tsv")
    file("${id}_daterandom.png")

    """
    nanopath phybeast utils date-random-test --dates $meta_data --alignment $alignment --tree $tree --replicates $params.replicates
    nanopath phybeast utils date-random-test --rate_file rates.tsv --clock_rate_file $rate
    mv date_random_test/date_randomisation_test.png ${id}_daterandom.png
    mv rates.tsv ${id}_rates.tsv
    """

}

process VariantSites {

    label "phybeast"
    tag { "$tree" }

    publishDir "${params.outdir}/alignment", mode: "copy"

    input:
    tuple val(id), file(alignment)

    output:
    tuple val(id), file("${id}.noref.variants.fasta")

    """
    nanopath phybeast utils remove-reference -a $alignment -o aln.noref.fasta
    nanopath phybeast utils remove-invariant -a aln.noref.fasta -o ${id}.noref.variants.fasta
    """

}