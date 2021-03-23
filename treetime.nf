process TreeTime {

    label "treetime"

    publishDir "${params.outdir}/treetime", mode: "copy"

    input:
    tuple val(id), file(tree), file(alignment)
    file(meta_data)

    output:
    tuple val(id), file(tree), file(alignment), file("${id}_rate.txt")
    file("${id}_tt")


    """
    nanopath utils prepare-metadata -m $meta_data -p treetime -o treetime.meta
    treetime --tree $tree --aln $alignment --dates treetime.meta --branch-length-mode auto \
    --covariation --coalescent skyline --confidence --outdir ${id}_tt
    nanopath utils extract-rate -f ${id}_tt/molecular_clock.txt -p treetime -o ${id}_rate.txt
    """

}