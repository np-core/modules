    process SnippyFastq {

        label "snippy"
        tag { id }

        input:
        tuple val(id), file(forward), file(reverse)
        file(reference)

        output:
        file("$id") // whole output folder

        """
        snippy --cpus $task.cpus --outdir $id --prefix $id --reference $reference --R1 $forward --R2 $reverse
        """

    }

    process SnippyFasta {

        label "snippy"
        tag { id }

        input:
        tuple val(id), file(fasta)
        file(reference)

        output:
        file("$id") // whole output folder

        """
        snippy --cpus $task.cpus --outdir $id --prefix $id --reference $reference --ctgs $fasta
        """

    }

    process SnippyCore {

        label "snippy"
        tag { "CoreAlignment" }

        publishDir "${params.outdir}/alignment", mode: "copy", pattern: "*.fasta"

        input:
        file(snippy_outputs)  // collected list of snippy output directories
        file(reference)

        output:
        file("wgs.core.fasta")

        """
        snippy-core --ref $reference --prefix core $snippy_outputs
        mv core.aln snp.core.fasta
        snippy-clean_full_aln core.full.aln > wgs.core.fasta
        """

    }