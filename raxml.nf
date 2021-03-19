process RAxML {

    label "raxml"
    tag { "$params.tree_model" }

    publishDir "${params.outdir}/phylogeny", mode: "copy"

    input:
    tuple val(id), file(alignment)

    output:
    tuple val(id), file("${id}.newick")

    """
    raxml-ng --msa $alignment --model $params.raxml_model $params.raxml_params --threads $task.cpus --prefix rax --force
    mv rax.raxml.bestTree ${id}.newick
    """

}
