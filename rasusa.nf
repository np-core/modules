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

process RasusaMultiTraining {
    
    tag { "${id} - ${model_name}: ${coverage}x" }
    label "rasusa"

    publishDir "$params.outdir/rasusa", mode: "copy"

    input:
    tuple val(model_name), val(id), file(fq), file(snippy_vcf)
    each coverage

    output:
    tuple val(model_name), val(id), file("${id}_${coverage}.fq"), file(snippy_vcf)
    
    """
    rasusa -c $coverage -g $params.genome_size -i $fq > ${id}_${coverage}.fq
    """
}