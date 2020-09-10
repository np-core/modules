
process Racon {

    tag { id }
    label "racon"

    input:
    tuple val(id), file(assembly), file(fastq)

    output:
    tuple val(id), file("${id}.racon.fasta"), file(fastq)

    script:
    """
    minimap2 -x map-ont -t $task.cpus $assembly $fastq > assembly_1.paf
    racon -m 8 -x -6 -g -8 -w 500 -t $task.cpus $fastq assembly_1.paf $assembly > assembly_consensus_1.fasta
    minimap2 -x map-ont -t $task.cpus assembly_consensus_1.fasta $fastq > assembly_2.paf
    racon -m 8 -x -6 -g -8 -w 500 -t $task.cpus $fastq assembly_2.paf assembly_consensus_1.fasta > assembly_consensus_2.fasta
    minimap2 -x map-ont -t $task.cpus assembly_consensus_2.fasta $fastq > assembly_3.paf
    racon -m 8 -x -6 -g -8 -w 500 -t $task.cpus $fastq assembly_3.paf assembly_consensus_2.fasta > assembly_consensus_3.fasta
    minimap2 -x map-ont -t $task.cpus assembly_consensus_3.fasta $fastq > assembly_4.paf
    racon -m 8 -x -6 -g -8 -w 500 -t $task.cpus $fastq assembly_4.paf assembly_consensus_3.fasta > ${id}.racon.fasta
    """
}