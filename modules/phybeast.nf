process DateRandomisation {

    label "phybeast"
    tag { "$tree" }

    publishDir "${params.outdir}", mode: "copy"

    input:
    file(tree)
    file(rate)
    file(meta_data)
    file(alignment)

    output:
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

    publishDir "${params.outdir}", mode: "copy"

    input:
    file(alignment)

    output:
    file("snps.noref.variants.fasta")

    """
    nanopath phybeast utils remove-reference -a $alignment -o snps.noref.fasta
    nanopath phybeast utils remove-invariant -a snps.noref.fasta -o snps.noref.variants.fasta
    """

}