process Beast {

    // BEAST2
    
    label "$beast_label"
    tag { "BEAST2" }

    publishDir "${params.outdir}/beast", mode: "copy"

    input:
    tuple val(id), file(xml)
    val(beagle_params)
    val(beast_label)

    output:
    tuple val(id), file("${id}.*")

    """
    beast -threads $task.cpus ${beagle_params} ${params.beast_params} $xml
    """

}
