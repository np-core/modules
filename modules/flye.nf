process Flye {
    
    tag { id }
    label "flye"

    publishDir "$params.outdir/ont/assembly", mode: "copy", pattern: "*.fasta"
    publishDir "$params.outdir/ont/assembly", mode: "copy", pattern: "*.gfa"
    publishDir "$params.outdir/ont/assembly", mode: "copy", pattern: "*.txt"

    input:
    tuple val(id), file(fq)

    output:
    tuple val(id), file("${id}.fasta")

    """
    flye --nano-raw $fq $params.assembly_options -t $task.cpus -o assembly
    mv assembly/assembly_info.txt ${id}.txt
    mv assembly/assembly.fasta ${id}.fasta
    mv assembly/assembly_graph.gfa ${id}.gfa
    """

}

process MetaFlye {
    
    tag { "$id" }
    label "flye"

    publishDir "$params.outdir/assembly/metaflye", mode: "copy", pattern: "*.fasta"
    publishDir "$params.outdir/assembly/metaflye", mode: "copy", pattern: "*.gfa"
    publishDir "$params.outdir/assembly/metaflye", mode: "copy", pattern: "*.txt"

    input:
    tuple val(id), file(fq)

    output:
    tuple val(id), file("${id}.fasta") optional true  // can fail to produce meaningful assembly on small files, fail ok
    tuple file("${id}.gfa"), file("${id}.txt") optional true

    """
    flye --nano-raw $fq --meta $params.assembly_options -t $task.cpus -o assembly
    
    if [ -f assembly/assembly.fasta ]; then
        mv assembly/assembly_info.txt ${id}.txt
        mv assembly/assembly.fasta ${id}.fasta
        mv assembly/assembly_graph.gfa ${id}.gfa
    fi
    """

}