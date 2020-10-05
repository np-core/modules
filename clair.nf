process ClairVariants {

    label "clair"
    tag { "$id" }
    
    memory { 8.GB * task.attempt }

    errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.outdir}/clair", mode: "copy", pattern: "${id}.vcf"
    publishDir "${params.outdir}/clair", mode: "copy", pattern: "${id}.txt"

    input:
    tuple val(id), file(bam), file(bai)
    file(reference) 

    output:
    tuple val(id), file("${id}.vcf"), file("${id}.txt")
    tuple val(id), file(bam), file(bai)

    """
    samtools faidx $reference
    np phybeast utils print-header --fasta $reference | while read -r contig ; do
    echo "Processing contig: \$contig"
    clair callVarBam --chkpnt_fn ${params.clair_model} \
                     --ref_fn $reference \
                     --bam_fn $bam \
                     --sampleName $id \
                     --minCoverage 1 \
                     --threads $task.cpus \
                     --call_fn ${id}.\${contig}.clair.vcf \
                     --ctgName \$contig \
                     ${params.clair_haploid}
    done

    vcfcat ${id}.*.clair.vcf | bcftools sort -m 2G -o ${id}.vcf 

    pysamstats -t variation_strand $bam -f $reference > ${id}.txt
    
    """

}

process ClairVariantsTraining {

    label "clair"
    tag { "$id" }
    
    memory { 8.GB * task.attempt }

    errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.outdir}/clair/${model_name}/${reference.baseName}", mode: "copy", pattern: "*_${coverage}.vcf"
    publishDir "${params.outdir}/clair/${model_name}/${reference.baseName}", mode: "copy", pattern: "*_${coverage}.txt"

    input:
    tuple val(model_name), val(coverage), val(ids), file(reference), file("bams/*"), file("bams/*")

    output:
    tuple val(model_name), file("*.vcf"), file("*.txt")


    """
    samtools faidx $reference
    
    for i in $ids; do
        np phybeast utils print-header --fasta $reference | while read -r contig ; do
        echo "Processing contig: \$contig"
        clair callVarBam --chkpnt_fn ${params.clair_model} \
                        --ref_fn $reference \
                        --bam_fn bams/\${i}_${coverage}.bam \
                        --sampleName \${i}_${coverage} \
                        --minCoverage 1 \
                        --threads $task.cpus \
                        --call_fn \${i}_${coverage}.\${contig}.clair.vcf \
                        --ctgName \$contig \
                        ${params.clair_haploid}
        done

        vcfcat \${i}_${coverage}.*.clair.vcf | bcftools sort -m 2G -o \${i}_${coverage}.vcf 

        pysamstats -t variation_strand bams/\${i}_${coverage}.bam -f $reference > \${i}_${coverage}.txt               
    done
    """

}