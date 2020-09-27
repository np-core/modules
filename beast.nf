process BeastCPU {

    // BEAST2 CPU
    
    label "beast"
    tag { "BEAST2" }

    publishDir "${params.outdir}/beast", mode: "copy"

    input:
    tuple val(id), file(xml)
    val(beagle_params)

    output:
    tuple val(id), file("${id}.*")

    """
    beast -threads $task.cpus ${beagle_params} ${params.beast_params} $xml
    """

}

process BeastGPU {

    // BEAST2 GPU
    
    label "beast_gpu"
    tag { "BEAST2" }

    publishDir "${params.outdir}/beast", mode: "copy"

    input:
    tuple val(id), file(xml)
    val(beagle_params)

    output:
    tuple val(id), file("${id}.*")

    """
    beast -threads $task.cpus ${beagle_params} ${params.beast_params} $xml
    """

}

