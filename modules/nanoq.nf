process Nanoq {
    
    // Nanoq read stream filter and summary statistics

    tag { id }
    label "nanoq"

    publishDir "$params.outdir/ont/qc", mode: "copy", pattern: "*.txt"

    input:
    tuple val(id), file(fq)

    output:
    tuple val(id), file("${id}.filtered.fq")

    """
    if [[ $fq == *.gz ]]
    then
        zcat $fq | nanoq -l $params.length -q $params.quality > ${id}.filtered.fq 2> ${id}.filtered.stats.txt
    else
        nanoq -f $fq -l $params.length -q $params.quality > ${id}.filtered.fq 2> ${id}.filtered.stats.txt
    fi
    """

}

process NanoqStatistics {
    
    // Fast variant to compuute statistics only

    tag { id }
    label "nanoq"

    publishDir "$params.outdir/ont/qc", mode: "copy", pattern: "*.txt"

    input:
    tuple val(id), file(fq)

    output:
    tuple val(id), file("${id}.stats.txt")

    """
    if [[ $fq == *.gz ]]
    then
        zcat $fq | nanoq -l $params.length -q $params.quality > /dev/null 2> ${id}.stats.txt
    else
        nanoq -f $fq -l $params.length -q $params.quality > /dev/null 2> ${id}.stats.txt
    fi
    """

}

// Pathogen variants

process NanoqOnline {
    
    // Nanoq read stream filter and summary statistics batch-wise

    tag { "Batch $batch - $id" }
    label "nanoq"

    publishDir "$params.outdir/fastq/nanoq", mode: "copy", pattern: "${id}.${batch}.filtered.stats.txt"
    publishDir "$params.outdir/fastq/nanoq", mode: "symlink", pattern: "${id}.${batch}.filtered.fq"

    input:
    tuple val(id), file(fq), val(batch)

    output:
    tuple val(id), file("${id}.${batch}.filtered.fq"), val(batch)
    file("${id}.${batch}.filtered.stats.txt")

    """
    if [[ $fq == *.gz ]]
    then
        zcat $fq | nanoq -l $params.min_length -q $params.min_quality > ${id}.${batch}.filtered.fq 2> ${id}.${batch}.filtered.stats.txt
    else
        nanoq -f $fq -l $params.min_length -q $params.min_quality > ${id}.${batch}.filtered.fq 2> ${id}.${batch}.filtered.stats.txt
    fi
    """

}

process NanoqOffline {
    
    // Nanoq read stream filter and summary statistics

    tag { id }
    label "nanoq"

    publishDir "$params.outdir/fastq/nanoq", mode: "symlink", pattern: "${id}.filtered.fq"
    publishDir "$params.outdir/fastq/nanoq", mode: "copy", pattern: "${id}.filtered.stats.txt"


    input:
    tuple val(id), file(fq)

    output:
    tuple val(id), file("${id}.filtered.fq")
    file("${id}.filtered.stats.txt")


    """
    if [[ $fq == *.gz ]]
    then
        zcat $fq | nanoq -l $params.min_length -q $params.min_quality > ${id}.filtered.fq 2> ${id}.filtered.stats.txt
    else
        nanoq -f $fq -l $params.min_length -q $params.min_quality > ${id}.filtered.fq 2> ${id}.filtered.stats.txt
    fi
    """

}