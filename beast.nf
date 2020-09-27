process Beast {

    // BEAST2
    
    label "beast2"
    tag { "" }

    publishDir "${params.outdir}/beast", mode: "copy"

    input:
    tuple val(id), file(xml)
    val(beagle_params)

    output:
    tuple val(id), file("${id}.*")

    """
    beast -threads $task.cpus ${beagle_params} ${params.beagle_order} ${params.beast_params} $xml
    """

}