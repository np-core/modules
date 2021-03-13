process Shovill {

    tag { id }
    label "shovill"

    publishDir "$params.outdir/illumina/assembly", mode: "copy", pattern: "*.fasta"
    publishDir "$params.outdir/illumina/assembly", mode: "copy", pattern: "*.gfa"

    input:
    tuple val(id), file(forward), file(reverse)

    output:
    tuple val(id), file("${id}.assembly.fasta")

    """
    shovill --R1 $forward --R2 $reverse --cpus $task.cpus --ram $task.memory \
    --depth $params.depth --assembler $params.assembler --outdir $id --force
    mv ${id}/contigs.fa ${id}.assembly.fasta
    """

}