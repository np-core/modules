process Rasusa {
    
    tag { id }
    label "ont"

    publishDir "$params.outdir/rasusa", mode: "copy"

    input:
    tuple val(id), file(fq)

    output:
    tuple val(id), file("${id}.rasusa.fq")
    
    script:

    if ( params.subsample > 0 )
        
        """
        rasusa -c $params.subsample -g $params.genome_size -i $fq > ${id}.rasusa.fq
        """
    else 
        """
        cp $fq ${id}.rasusa.fq
        """

}

process RasusaMulti {
    
    tag { id }
    label "ont"

    publishDir "$params.outdir/rasusa", mode: "copy"

    input:
    tuple val(id), file(fq)
    each coverage

    output:
    tuple val(id), file("${id}_${coverage}.fq")
    
    """
    rasusa -c $coverage -g $params.genome_size -i $fq > ${id}_${coverage}.fq
    """
}