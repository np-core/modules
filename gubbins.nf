process Gubbins {

    label "gubbins"
    tag { "Gubbins" }

    publishDir "${params.outdir}", mode: "copy", pattern: "snp.core.gubbins.fasta"

    input:
    file(clean_alignment)

    output:
    file("snp.core.gubbins.fasta")

    """
    run_gubbins.py -p gubbins --threads $task.cpus $clean_alignment $params.gubbins_params
    snp-sites -c gubbins.filtered_polymorphic_sites.fasta > snp.core.gubbins.fasta
    """

}