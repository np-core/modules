process RAxML {

    label "raxml"
    tag { "$params.tree_model" }

    publishDir "${params.outdir}/phylogeny", mode: "copy"

    input:
    file(alignment)

    output:
    file("tree.newick")

    """
    raxml-ng --msa $alignment --model $params.raxml_model $params.raxml_params --threads $task.cpus --prefix rax --force
    mv rax.raxml.bestTree tree.newick
    """

}
