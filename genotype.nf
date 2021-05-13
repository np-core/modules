 process AssemblyGenotype {

    tag { id }
    label "genotype"

    publishDir "$params.outdir/${params.tag}/genotypes", mode: "copy", pattern: "${id}.${params.tag}.tab"

    input:
    tuple val(id), file(assembly)

    output:
    file("${id}.${params.tag}.tab")

    script:

    if (params.kpneumoniae)
        
        """
        kleborate -a $assembly --all -o ${id}.${params.tag}.tab
        """

    else if (params.saureus)
        
        """
        sccion type -a $assembly > ${id}.${params.tag}.tab
        """
    
    else

        """
        mlst $assembly >> ${id}.${params.tag}.tab
        abricate --db vfdb $assembly >> ${id}.${params.tag}.tab
        """

}

 process Mykrobe {

    tag { id }
    label "genotype"

    publishDir "$params.outdir/${params.tag}/genotypes", mode: "copy", pattern: "${id}.json"

    input:
    tuple val(id), file(forward), file(reverse)

    output:
    file("${id}.json")

    script:

    """
    mykrobe predict --sample $id --species $params.mykrobe_species $params.mykrobe_other --out ${id}.json --format json -i $forward $reverse
    """

}