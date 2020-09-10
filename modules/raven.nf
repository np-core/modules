process RavenRegexOnline  {
  
    tag { "Batch: $batch  - MAGs" }
    label "raven"

    publishDir "$params.outdir/assembly/raven", mode: "copy", pattern: "${id}.${batch}.*.fa"

    input:
    tuple val(id), file(fq), val(batch)
    file(assembly_update)  // this exists if read files were linked beforehand - LinkAssemblyOnline

    output:
    tuple val(id), file("${id}.${batch}.*.fa"), val(batch)  // pipe this into flat map to trigger kraken for each regex group mag assembly (barcode)
    
    when:
    batch % params.assembly_every == 0  // trigger this assembly only every

    shell:

    // TODO: when Raven fail (no files, regex groups)

    """
    np pathogen raven-mag -p $assembly_update --file_glob "*" --regex "$params.read_regex" --outdir workdir --raven_args "$params.assembly_options"
    
    for f in workdir/*.fa; do
        group=\$(basename \$f .fa)
        mv \$f ${id}.${batch}.\${group}.fa
    done
    """
}

process RavenRegex {
  
    tag { "$id" }
    label "raven"

    publishDir "$params.outdir/assembly/raven", mode: "copy", pattern: "${id}.*.fa"

    input:
    tuple val(id), file(fq)
    file(assembly_update)

    output:
    tuple val(id), file("${id}.*.fa")  // pipe this into flat map to trigger kraken for each regex group mag assembly (barcode)
    
    shell:

    """
    np pathogen raven-mag -p $assembly_update --file_glob "*" --regex "$params.read_regex" --outdir workdir --raven_args "$params.assembly_options"
    for f in workdir/*.fa; do
        group=\$(basename \$f .fa)
        mv \$f ${id}.\${group}.fa
    done
    """
}