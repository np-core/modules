process ClairVariants {

    label "clair"
    tag { "$id" }
    
    memory { 4.GB * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
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
                     --haploid_sensitive
    done

    vcfcat ${id}.*.clair.vcf | bcftools sort -m 2G -o ${id}.vcf 

    pysamstats -t variation_strand $bam -f $reference > ${id}.txt
    
    """

}