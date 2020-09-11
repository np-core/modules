
process Guppy {

    tag { id }
    label "guppy"

    publishDir "$params.outdir/guppy", pattern: "*.telemetry"
    publishDir "$params.outdir/guppy", pattern: "*.summary"

    input:
    tuple val(id), file(path)

    output:
    tuple val(id), file("${id}.fq")
    tuple file("${id}.summary"), file("${id}.telemetry")

    """
    mkdir fast5

    if [[ ($path == *.tar.gz) || ($path == *.tar) ]]
    then
        tar -xf $path -C fast5
    else
        mv $path fast5
    fi

    guppy_basecaller -i fast5 -s $id -c $params.guppy_model -x "$params.devices" $params.guppy_params \
        --gpu_runners_per_device $params.runners_per_device --chunk_size $params.chunk_size \
        --chunks_per_runner $params.chunks_per_runner --num_callers $params.num_callers -r

    cat $id/*.fastq > ${id}.fq

    mv $id/sequencing_telemetry.js ${id}.telemetry
    mv $id/sequencing_summary.txt ${id}.summary
    """

}