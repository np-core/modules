process TreeTime {

    label "treetime"

    publishDir "${params.outdir}", mode: "copy"

    input:
    file(tree)
    file(meta_data)
    file(alignment)

    output:
    file("rate.txt")
    file("treetime")


    """
    nanopath phybeast utils prepare-metadata -m $meta_data -p treetime -o treetime.meta
    treetime --tree $tree --aln $alignment --dates treetime.meta --branch-length-mode auto \
    --covariation --coalescent skyline --confidence --outdir treetime
    nanopath phybeast utils extract-rate -f treetime/molecular_clock.txt -p treetime -o rate.txt
    """

}