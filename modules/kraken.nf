process Kraken {

    tag { id }
    label "kraken"

    publishDir "$params.outdir/kraken/$db", mode: "copy", pattern: "*.kraken"
    publishDir "$params.outdir/kraken/$db", mode: "copy", pattern: "*.bracken"
    publishDir "$params.outdir/kraken/$db", mode: "copy", pattern: "*.report"
    publishDir "$params.outdir/fastq", pattern: "$fq"  // symlink the fastq (e.g. in online run 4000 reads default) for unorthodox online coverage assessment triggered after this process

    input:
    tuple val(id), file(fq)
    each file(db)

    output:
    tuple val(id), file("${id}.kraken"), file("${id}.kraken.report"), file("${id}.bracken"), file("${id}.bracken.report")

    """
    kraken2 --db $db --threads $task.cpus --output ${id}.kraken --report ${id}.kraken.report $fq
    bracken -d $db -i ${id}.kraken.report -o ${id}.bracken -w ${id}.bracken.report -r $params.bracken_length -l $params.bracken_level -t $params.bracken_threshold
    """

}

process KrakenOnline {

    tag { "Batch $batch - $db" }
    label "kraken"

    publishDir "$params.outdir/kraken/$db", mode: "copy", pattern: "*.kraken"
    publishDir "$params.outdir/kraken/$db", mode: "copy", pattern: "*.bracken"
    publishDir "$params.outdir/kraken/$db", mode: "copy", pattern: "*.report"

    input:
    tuple val(id), file(fq), val(batch)
    each file(db)

    output:
    tuple val(id), file("${id}.${batch}.kraken"), file("${id}.${batch}.kraken.report"), file("${id}.${batch}.bracken"), file("${id}.${batch}.bracken.report")

    """
    kraken2 --db $db --threads $task.cpus --output ${id}.${batch}.kraken --report ${id}.${batch}.kraken.report $fq
    bracken -d $db -i ${id}.${batch}.kraken.report -o ${id}.${batch}.bracken -w ${id}.${batch}.bracken.report -r $params.bracken_length -l $params.bracken_level -t $params.bracken_threshold
    """

}

process KrakenAssemblyOnline {

    tag { "Batch $batch - $db" }
    label "kraken"

    publishDir "$params.outdir/assembly/kraken/$db", mode: "copy", pattern: "*.kraken"
    publishDir "$params.outdir/assembly/kraken/$db", mode: "copy", pattern: "*.bracken"
    publishDir "$params.outdir/assembly/kraken/$db", mode: "copy", pattern: "*.report"

    input:
    tuple val(id), file(fq), val(batch)
    each file(db)

    output:
    tuple val(id), file("${id}.${batch}.kraken"), file("${id}.${batch}.kraken.report"), file("${id}.${batch}.bracken"), file("${id}.${batch}.bracken.report")

    """
    kraken2 --db $db --threads $task.cpus --output ${id}.${batch}.kraken --report ${id}.${batch}.kraken.report $fq
    bracken -d $db -i ${id}.${batch}.kraken.report -o ${id}.${batch}.bracken -w ${id}.${batch}.bracken.report -r $params.bracken_length -l $params.bracken_level -t $params.bracken_threshold
    """

}


process KrakenIllumina {

    tag { id }
    label "kraken"

    publishDir "$params.outdir/kraken/$db", mode: "copy", pattern: "*.kraken"
    publishDir "$params.outdir/kraken/$db", mode: "copy", pattern: "*.bracken"
    publishDir "$params.outdir/kraken/$db", mode: "copy", pattern: "*.report"

    input:
    tuple val(id), file(forward), file(reverse)
    each file(db)

    output:
    tuple val(id), file("${id}.kraken"), file("${id}.kraken.report"), file("${id}.bracken"), file("${id}.bracken.report")

    """
    kraken2 --db $db --threads $task.cpus --output ${id}.kraken --report ${id}.kraken.report --paired $forward $reverse
    bracken -d $db -i ${id}.kraken.report -o ${id}.bracken -w ${id}.bracken.report -r $params.bracken_length -l $params.bracken_level -t $params.bracken_threshold
    """

}

