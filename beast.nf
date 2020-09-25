process Beast {

    // BEAST2
    
    label "beast2"
    tag { "" }

    publishDir "${params.outdir}/beast", mode: "copy"

    input:
    tuple val(id), file(xml)
    val(beagle_params)

    output:
    tuple id, file("${id}.*")

    """
    beast ${beagle_params} ${params.beagle_order} ${params.beat_params} $xml
    """

}