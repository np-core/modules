process MegalodonVariants {

    label "megalodon"
    tag { "Megalodon: $panel - $barcode" }

    publishDir "${params.outdir}/megalodon", mode: "copy"

    input:
    file(reference)
    tuple val(port), val(panel), val(barcode), file(path)

    output:
    file("${panel}_${barcode}")

    """
    megalodon --guppy-server-path $params.guppy_server_path \
            --output-directory ${panel}_${barcode} \
            --outputs variants \
            --reference $reference \
            --haploid \
            --variant-filename $params.candidates \
            --devices $params.devices \
            --processes $task.cpus \
            --guppy-server-port $port \
            --guppy-params "$params.guppy_params" \
            --guppy-config $params.guppy_config 
            $path 
    """

}