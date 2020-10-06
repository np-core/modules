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
    tag { "$model_name - $id - $reference" }
    
    memory { 8.GB * task.attempt }

    errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.outdir}/clair/${model_name}/${reference.baseName}", mode: "copy", pattern: "${id}_${coverage}.vcf"
    publishDir "${params.outdir}/clair/${model_name}/${reference.baseName}", mode: "copy", pattern: "${id}_${coverage}.txt"

    input:
    tuple val(model_name), val(id), val(coverage), file(reference), file(bam), file(bai), file(snippy_vcf)

    output:
    tuple val(model_name), val("${reference.baseName}"), file("${id}_${coverage}.vcf"), file("${id}_${coverage}.txt"), file(snippy_vcf)


    """
    samtools faidx $reference
    
    np phybeast utils print-header --fasta $reference | while read -r contig ; do
    echo "Processing contig: \$contig"
    clair callVarBam --chkpnt_fn ${params.clair_model} \
                    --ref_fn $reference \
                    --bam_fn $bam \
                    --sampleName ${id}_${coverage} \
                    --minCoverage 1 \
                    --threads $task.cpus \
                    --call_fn ${id}_${coverage}.\${contig}.clair.vcf \
                    --ctgName \$contig \
                    ${params.clair_haploid}
    done

    vcfcat ${id}_${coverage}.*.clair.vcf | bcftools sort -m 2G -o ${id}_${coverage}.vcf 

    pysamstats -t variation_strand $bam -f $reference > ${id}_${coverage}.txt               

    """

}