process MegalodonVariants {

    label "megalodon"
    tag { "Megalodon: $id" }

    publishDir "${params.outdir}/megalodon", mode: "copy"

    input:
    file(reference)
    file(candidates)
    tuple val(id), file(path)

    output:
    file("${id}")

    """
    mv $path megalodon_in
    megalodon --guppy-server-path $params.guppy_server_path \
            --output-directory ${id} \
            --outputs variants \
            --reference $reference \
            --haploid \
            --variant-filename $candidates \
            --devices $params.gpu_devices \
            --processes $task.cpus \
            --guppy-params "$params.guppy_params" \
            --guppy-config "$params.guppy_config" \
            $params.megalodon_params megalodon_in
    """

}

process MegalodonVariantsPanels {

    label "megalodon"
    tag { "Megalodon: $panel - $barcode" }

    publishDir "${params.outdir}/megalodon", mode: "copy"

    input:
    file(reference)
    file(candidates)
    tuple val(panel), val(barcode), file(path)

    output:
    file("${panel}_${barcode}")

    """
    megalodon --guppy-server-path $params.guppy_server_path \
            --output-directory ${panel}_${barcode} \
            --outputs variants \
            --reference $reference \
            --haploid \
            --variant-filename $candidates \
            --devices $params.devices \
            --processes $task.cpus \
            --guppy-params "$params.guppy_params" \
            --guppy-config $params.guppy_config \
            $params.megalodon_params $path
    """

}

