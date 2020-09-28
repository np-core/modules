process MinimapONT {

    label "minimap2"
    tag { "$id" }

    publishDir "${params.outdir}/minimap2", mode: "symlink"

    input:
    tuple val(id), file(fastq)
    file(reference)

    output:
    tuple val(id), file("${id}.bam"), file("${id}.bam.bai")

    """
    minimap2 -t $task.cpus -ax map-ont $reference $fastq | samtools sort | samtools view -Sb > ${id}.bam
    samtools index ${id}.bam
    """

}