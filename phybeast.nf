process DateRandomisation {

    label "phybeast"
    tag { "$tree" }

    publishDir "${params.outdir}/date_randomisation", mode: "copy"

    input:
    file(tree)
    file(rate)
    file(meta_data)
    file(alignment)

    output:
    file("rates.tsv")
    file("date_randomisation_test.png")

    """
    nanopath phybeast utils date-random-test --dates $meta_data --alignment $alignment --tree $tree --replicates $params.replicates
    nanopath phybeast utils date-random-test --rate_file rates.tsv --clock_rate_file $rate
    mv date_random_test/date_randomisation_test.png date_randomisation_test.png
    """

}

process VariantSites {

    label "phybeast"
    tag { "$tree" }

    publishDir "${params.outdir}/alignment", mode: "copy"

    input:
    tuple val(id), file(alignment)

    output:
    file("${id}.noref.variants.fasta")

    """
    nanopath phybeast utils remove-reference -a $alignment -o aln.noref.fasta
    nanopath phybeast utils remove-invariant -a aln.noref.fasta -o aln.noref.variants.fasta
    """

}