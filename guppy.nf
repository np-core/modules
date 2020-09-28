
process Guppy {

    tag { id }
    label "guppy"

    publishDir "$params.outdir/guppy", mode: 'copy', pattern: "*.telemetry"
    publishDir "$params.outdir/guppy", mode: 'copy', pattern: "*.summary"
    publishDir "$params.outdir/guppy", mode: 'copy', pattern: "*.fq"
    publishDir "$params.outdir/guppy", mode: 'symlink', pattern: "${id}"

    input:
    tuple val(id), file(path)

    output:
    tuple val(id), file("${id}.fq")
    tuple file("${id}.summary"), file("${id}.telemetry"), file("${id}")

    """

    if [[ -f $path ]]; then
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
        mv $path fast5_in
    else
        echo "Error in parsing input"
        exit 1
    fi

    guppy_basecaller -i fast5_in -s $id -c $params.guppy_model -d "$params.guppy_data" -x "$params.guppy_devices" $params.guppy_params \
        --gpu_runners_per_device $params.runners_per_device --chunk_size $params.chunk_size \
        --chunks_per_runner $params.chunks_per_runner --num_callers $params.num_callers \
        --cpu_threads_per_caller $params.cpu_threads_per_caller

    cat $id/*.fastq > ${id}.fq

    mv $id/sequencing_telemetry.js ${id}.telemetry
    mv $id/sequencing_summary.txt ${id}.summary
    """

}


process GuppyBatch {

    tag { id }
    label "guppy"

    publishDir "$params.outdir/guppy", mode: 'copy', pattern: "*.telemetry"
    publishDir "$params.outdir/guppy", mode: 'copy', pattern: "*.summary"
    publishDir "$params.outdir/guppy", mode: 'copy', pattern: "*.fq"
    publishDir "$params.outdir/guppy", mode: 'symlink', pattern: "${id}"

    input:
    tuple val(id), file(path)

    output:
    tuple val(id), file("${id}.fq")
    tuple file("${id}.summary"), file("${id}.telemetry"), file("${id}")

    """
    mkdir fast5_in
    mv $path fast5_in

    guppy_basecaller -i fast5_in -s $id -c $params.guppy_model -d "$params.guppy_data" -x "$params.guppy_devices" $params.guppy_params \
        --gpu_runners_per_device $params.runners_per_device --chunk_size $params.chunk_size \
        --chunks_per_runner $params.chunks_per_runner --num_callers $params.num_callers \
        --cpu_threads_per_caller $params.cpu_threads_per_caller

    cat $id/*.fastq > ${id}.fq

    mv $id/sequencing_telemetry.js ${id}.telemetry
    mv $id/sequencing_summary.txt ${id}.summary
    """

}