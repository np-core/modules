
process Bonito {

    tag { id }
    label "bonito"

    publishDir "$params.outdir/bonito", mode: 'copy', pattern: "*.summary"
    publishDir "$params.outdir/bonito", mode: 'copy', pattern: "*.fq"

    input:
    tuple val(id), file(path)

    output:
    tuple val(id), file("${id}.fq")
    file("${id}.summary")

    """

    if [[ -f $path ]]; then
        # if the path variable is a single file
        if [[ ($path == *.tar.gz) || ($path == *.tar) ]]; then
            # if it is archived
            mkdir fast5_in
            tar -xf $path -C fast5_in
        else
            # if it is a single file
            mkdir fast5_in
            mv $path fast5_in
        fi
    elif [[ -d $path ]]; then
        # if the path variable is a single directory
        mkdir fast5_in
        mv $path fast5_in
    else
        echo "Error in parsing input"
        exit 1
    fi

    bonito basecall --device $params.bonito_device --fastq $params.bonito_model fast5_in > ${id}.fq
    mv summary.tsv ${id}.summary
    """

}


process BonitoBatch {

    tag { id }
    label "bonito"

    publishDir "$params.outdir/guppy", mode: 'copy', pattern: "*.summary"
    publishDir "$params.outdir/guppy", mode: 'copy', pattern: "*.fq"

    input:
    tuple val(id), file(path)

    output:
    tuple val(id), file("${id}.fq")
    file("${id}.summary")

    """
    mkdir fast5_in
    mv $path fast5_in

    bonito basecall --device $params.bonito_device --fastq $params.bonito_model fast5_in > ${id}.fq
    mv summary.tsv ${id}.summary
    """

}