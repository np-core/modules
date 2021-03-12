process UnicyclerHybrid {

    label "unicycler"
    tag { "$id" }

    memory { params.forest_evaluate_mem * task.attempt }

    errorStrategy { task.exitStatus in 137..143 ? 'retry' : 'ignore' }
    maxRetries 5

    publishDir "${params.outdir}/hybrid/unicycler", mode: "copy"

    input:
    tuple val(id), file(fwd), file(rev), file(fq)

    output:    
    tuple val(id), file("${id}.fasta")

    """
    unicycler -1 $fwd -2 $rev -l $fq -o unicycler --threads $task.cpus --no_correct
    mv unicycler/assembly.fasta ${id}.fasta
    """

}

