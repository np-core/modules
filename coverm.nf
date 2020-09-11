process CoverM {

    tag { id }
    label "coverm"

    publishDir "$params.outdir/ont/qc", mode: "copy", pattern: "*.txt"

    input:
    tuple val(id), file(fastq)
    file(reference)
    
    output:
    tuple val(id), file("${id}.coverage.txt")

    """
    coverm genome -p $params.preset -r $reference -t $task.cpus --single $fastq coverm genome --single-genome -m mean > ${id}.coverage.txt
    """

}

process IlluminaCoverM {

    tag { id }
    label "coverm"

    publishDir "$params.outdir/ont/qc", mode: "copy", pattern: "*.txt"

    input:
    tuple val(id), file(fastq)
    file(reference)
    
    output:
    tuple val(id), file("${id}.coverage.txt")

    """
    coverm genome -p $params.preset -r $reference -t $task.cpus --single $fastq coverm genome --single-genome -m mean > ${id}.coverage.txt
    """

}

process OnlineCoverM {
    
    tag { id }
    label "coverm"

    publishDir "$params.outdir/ont/qc", mode: "copy", pattern: "*.txt"

    input:
    tuple val(id), file(fastq)
    file(reference)
    
    output:
    tuple val(id), file("${id}.coverage.txt") optional true

    // first group symlinked fastq files from kraken process in params.outdir/fastq by provided regex (usually: barcode\d\d)
    // and map against the ref input to get coverage for this selected reference, triggered after each kraken process and read 
    // file (usually: 400 reads in online run) is classified and references have been selected in feeding process

    """
    
    coverm genome -p $params.preset -r $reference -t $task.cpus --single $fastq coverm genome --single-genome -m mean > ${id}.coverage.txt
    """

}