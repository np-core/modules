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
    file("${id}_unicycler")

    """
    unicycler -1 $fwd -2 $rev -l $fq -o ${id}_unicycler --threads $task.cpus --no_correct
    """

}

