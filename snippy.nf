    process SnippyFastq {

        label "snippy"
        tag { id }

        publishDir "${params.outdir}/snippy", mode: "symlink", pattern: "$id"

        input:
        tuple val(id), file(forward), file(reverse)
        file(reference)

        output:
        file("$id") // whole output folder

        """
        snippy --cpus $task.cpus --outdir $id --prefix $id --reference $reference --R1 $forward --R2 $reverse $params.snippy_params
        """

    }

    process SnippyFasta {

        label "snippy"
        tag { id }

        publishDir "${params.outdir}/snippy", mode: "symlink", pattern: "$id"

        input:
        tuple val(id), file(fasta)
        file(reference)

        output:
        file("$id") // whole output folder

        """
        snippy --cpus $task.cpus --outdir $id --prefix $id --reference $reference --ctgs $fasta $params.snippy_params
        """

    }

    process SnippyCore {

        label "snippy"
        tag { "SnippyCore" }

        publishDir "${params.outdir}/snippy_core", mode: "copy", pattern: "snp.core.vcf"
        publishDir "${params.outdir}/snippy_core", mode: "copy", pattern: "snp.core.fasta"
        publishDir "${params.outdir}/snippy_core", mode: "symlink", pattern: "wgs.core.fasta"

        input:
        file(snippy_outputs)  // collected list of snippy output directories
        file(reference)

        output:
        file("wgs.core.fasta")
        file("snp.core.fasta")
        file("snp.core.vcf")

        """
        snippy-core --ref $reference --prefix core $snippy_outputs
        mv core.aln snp.core.fasta
        mv core.vcf snp.core.vcf
        snippy-clean_full_aln core.full.aln > wgs.core.fasta
        """

    }


    process TrainingReferenceSnippy {
        
        // Calls the reference VCF from Illumina PE reads (reference) in training collections

        label "snippy"
        tag { id }

        publishDir "${params.outdir}/${ref}/polishers/snippy", mode: "symlink", pattern: "${id}_snippy/${id}.vcf"

        input:
        tuple val(model), val(id), file(ont), file(forward), file(reverse)
        each file(reference)

        output:
        tuple val(model), val(id), val(ref), file(reference), file(ont), file("${id}_snippy/${id}.vcf") 

        script:

        ref = reference.simpleName

        """
        snippy --cpus $task.cpus --outdir ${id}_snippy --prefix $id --reference $reference --R1 $forward --R2 $reverse $params.snippy_params
        """

    }