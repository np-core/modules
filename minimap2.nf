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
    tuple val(model_name), val(coverage), val(ids), file("ont_subset/*")
    each file(reference)

    output:
    tuple val(model_name), val(coverage), val(ids), file(reference), file("*.bam"), file("*.bam.bai")

    """
    for i in $ids; do
        echo \${i}_${coverage}
        minimap2 -t $task.cpus -ax map-ont $reference ont_subset/\${i}_${coverage}.fq | samtools sort | samtools view -Sb > \${i}_${coverage}.bam
        samtools index \${i}_${coverage}.bam
    done
    """

}