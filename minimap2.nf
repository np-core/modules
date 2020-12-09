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

process MinimapTraining {

    label "minimap2"
    tag { "$model_name - $id - $reference" }

    publishDir "${params.outdir}/${ref}/polishers/alignments", mode: "symlink"
 
    input:
    tuple val(model_name), val(id), val(ref), val(coverage), file(reference), file(fq_cov), file(snippy_vcf)

    output:
    tuple val(model_name), val(id), val(ref), val(coverage), file(reference), file("${id}_${coverage}.bam"), file("${id}_${coverage}.bam.bai"), file(snippy_vcf)

    """
    minimap2 -t $task.cpus -ax map-ont $reference $fq_cov | samtools sort | samtools view -Sb > ${id}_${coverage}.bam
    samtools index ${id}_${coverage}.bam
    """

}

process MinimapEvaluation {

    label "minimap2"
    tag { "$id - $reference" }

    memory { params.minimap_mem * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    publishDir "${params.outdir}/${reference.simpleName}/evaluation/${eval_set}", mode: "symlink"
 
    input:
    tuple val(eval_set), val(id), file(fq)
    each file(reference)

    output:
    tuple val(eval_set), val(id), file(reference), file("${id}.bam"), file("${id}.bam.bai")

    // Sometimes from Guppy basecalling there may be a mapping error of some reads (?)
    // Do not pipe as mapping will continue

    """
    minimap2 -t $task.cpus -ax map-ont $reference $fq > tmp.sam 
    samtools sort tmp.sam | samtools view -Sb > ${id}.bam
    samtools index ${id}.bam
    """

}