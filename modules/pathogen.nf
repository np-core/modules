process SelectReference {

    label "pathogen"
    tag { "#" }

    input:
    file(bracken_report)

    output:
    file("ref_*.fa") optional true 

    // Merges Bracken reports by Regex as specified (barcode or sample basename for example)


    // Collection must be setup in container from np pathogen setup --collection refseq --path /dbs/refseq [refseq, nctc]
    // Can fail to retrieve reference genome, therefore stops and does not process coverage or assembly data, 
    // needs to be checked in pathogen parser
    
    """
    nanopath pathogen utils select-reference --collections /dbs/refseq --bracken $params.outdir/ --coverage_threshold $params.coverage_threshold 
    """

}

process CoverageWindowONT {

    label "pathogen"
    tag { "#" }

    input:
    tuple val(id), file(fq)
    file(ref)

    output:
    file("${id}") optional true 

    
    """
    minimap2 -ax map-ont  $ref $fq | samtools sort -o aln.bam
    mosdepth -n --fast-mode --by $params.window_size-t $task.cpus --prefix $id aln.bam
    """

}

process CoverageWindowIllumina {

    label "pathogen"
    tag { "#" }



    input:
    tuple val(id), file(forward), file(reverse)
    file(ref)

    output:
    file("${id}") optional true 

    
    """
    minimap2 -ax sr $ref $forward $reverse | samtools sort -o aln.bam
    mosdepth -n --fast-mode --by $params.window_size $params.coverage_settings -t $task.cpus --prefix $id aln.bam
    """

}

process LinkAssemblyOnline {

    label "pathogen"
    tag { "Link reads to $params.assembly_update" }

    publishDir "$params.assembly_update", mode: "symlink"

    input:
    tuple val(id), file(fq), val(batch)

    output:
    tuple val(id), file(fq), val(batch)
    
    // All this does is symlink the fastq files so that the process that
    // derives from this has them available in the output dir

    """
    """

}

process FileGroupRegex {

    label "pathogen"
    tag { "Grouping files by regex: $params.read_regex"  }
    echo true

    input:
    val(reads) // from collect
    
    output:
    tuple val(id), file(fq)

    script:
    
    """
    echo $reads

    for i in $reads; do echo \$i; done
    """

}

process DeepHostRemovalOffline {
    
    // Use the global human diversity database 

    tag { "DeepHost: $params.host_database" }
    label "kraken2"

    publishDir "$params.outdir/fastq/host_removed", mode: "copy", pattern: "${id}.host.report"
    publishDir "$params.outdir/fastq/host_removed", mode: "copy", pattern: "${id}.host.kraken"
    publishDir "$params.outdir/fastq/host_removed", mode: "symlink", pattern: "${id}.nohost.fq"

    input:
    tuple val(id), file(fq)
    file(host_database)
    file(host_reference)

    output:
    tuple val(id), file("${id}.nohost.fq")
    tuple file("${id}.host.kraken"), file("${id}.host.report")

    shell:

    """
    kraken2 --db $params.host_database --report ${id}.host.report --unclassified-out ${id}.nohost.fq --output ${id}.host.kraken $fq
    """
}

process DeepHostRemovalOnline  {
    
    // Use the global human diversity database 

    tag { "Batch $batch - $id" }
    label "kraken2"

    publishDir "$params.outdir/fastq/host_removed", mode: "copy", pattern: "${id}.${batch}.host.report"
    publishDir "$params.outdir/fastq/host_removed", mode: "copy", pattern: "${id}.${batch}.host.kraken"
    publishDir "$params.outdir/fastq/host_removed", mode: "symlink", pattern: "${id}.${batch}.nohost.fq"

    input:
    tuple val(id), file(fq), val(batch)
    file(host_database)
    file(host_reference)

    output:
    tuple val(id), file("${id}.${batch}.nohost.fq"), val(batch)
    tuple file("${id}.${batch}.host.kraken"), file("${id}.${batch}.host.report")

    shell:

    """
    kraken2 --db $host_database --report ${id}.${batch}.host.report --unclassified-out ${id}.${batch}.nohost.fq --output ${id}.${batch}.host.kraken $fq
    """
}