process Pilon {

    tag { id }
    label "pilon"

    publishDir "$params.outdir/hybrid/assembly", mode: "copy", pattern: "${id}.*.fasta"

    input:
    tuple val(id), file(assembly), file(forward), file(reverse)

    output:
    tuple val(id), file("${id}.${params.tag}.hybrid.fasta")


    """
    minimap2 -ax sr $assembly $forward $reverse > aln.sam && \
        samtools view -S -b aln.sam > aln.bam && samtools sort aln.bam -o alignment1.bam && \
        samtools index alignment1.bam

    pilon --genome $assembly --frags alignment1.bam --outdir correction1 --changes
    
    minimap2 -ax sr correction1/pilon.fasta $forward $reverse > aln.sam && \
        samtools view -S -b aln.sam > aln.bam && samtools sort aln.bam -o alignment2.bam && \
        samtools index alignment2.bam

    pilon --genome correction1/pilon.fasta --frags alignment2.bam --outdir correction2 --changes && \
    mv correction2/pilon.fasta ${id}.${params.tag}.hybrid.fasta 
    """

}