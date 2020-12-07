process Rasusa {
    
    tag { id }
    label "rasusa"

    publishDir "$params.outdir/rasusa", mode: "copy"

    input:
    tuple val(id), file(fq)

    output:
    tuple val(id), file("${id}_rasusa.fq")
    
    """
    rasusa -c $params.coverage -g $params.genome_size -i $fq > ${id}_rasusa.fq
    """

}

process RasusaMulti {
    
    tag { "${id}: ${coverage}x" }
    label "rasusa"

    publishDir "$params.outdir/rasusa", mode: "copy"

    input:
    tuple val(id), file(fq)
    each coverage

    output:
    tuple val("${id}_${coverage}"), file("${id}_${coverage}.fq")
    
    """
    rasusa -c $coverage -g $params.genome_size -i $fq > ${id}_${coverage}.fq
    """
}

process RasusaTraining {
    
    tag { "${model_name} - ${id} - ${coverage}x" }
    label "rasusa"

    publishDir "$params.outdir/${ref}/polishers/subsets", mode: "copy", pattern: "${id}_${coverage}.fq"

    input:
    tuple val(model_name), val(id), val(ref), file(reference), file(fq), file(snippy_vcf)
    each coverage

    output:
    tuple val(model_name), val(id), val(ref), val(coverage), file(reference), file("${id}_${coverage}.fq"), file("${id}_${coverage}.ref.vcf")
    
    """
    rasusa -c $coverage -g $params.genome_size -i $fq > ${id}_${coverage}.fq
    cp $snippy_vcf ${id}_${coverage}.ref.vcf
    """
}