 process Genotype {

    tag { id }
    label "genotype"

    publishDir "$params.outdir/${params.tag}/genotypes", mode: "copy", pattern: "${id}.${params.tag}.tab"

    input:
    tuple val(id), file(assembly)

    output:
    file("${id}.${params.tag}.tab")

    when:
    params.kpneumoniae | params.saureus

    script:

    if (params.kpneumoniae)
        
        """
        kleborate -a $assembly --all -o ${id}.${params.tag}.tab
        """

    else if (params.saureus)
        
        """
        sccion type -a $assembly > ${id}.${params.tag}.tab
        """

}