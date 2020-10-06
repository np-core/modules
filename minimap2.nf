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

process MinimapMultiTraining {

    label "minimap2"
    tag { "$model_name" }

    publishDir "${params.outdir}/minimap2/${model_name}/${reference.baseName}", mode: "symlink"
 
    input:
    tuple val(model_name), val(id), val(coverage), file(fq_cov), file(snippy_vcf)
    each file(reference)

    output:
    tuple val(model_name), val(id), val(coverage), file(reference), file("${id}_${coverage}.bam"), file("${id}_${coverage}.bam.bai"), file(snippy_vcf)

    """
    minimap2 -t $task.cpus -ax map-ont $reference $fq_cov | samtools sort | samtools view -Sb > ${id}_${coverage}.bam
    samtools index ${id}_${coverage}.bam
    """

}