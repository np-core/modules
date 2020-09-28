process Minimap2ONT {

    label "minimap2"
    tag { "$id" }

    publishDir "${params.outdir}/medaka", mode: "copy"

    input:
    tuple val(id), file(fastq)
    file(reference)

    output:
    tuple val(id), file("${id}.bam"), file("${id}.bam.bai")

    """
    minimap2 -ax map-ont $reference $fastq | samtools sort | samtools view -Sb > ${id}.bam
    samtools index ${id}.bam
    """

}