    process SnippyFastq {

        label "snippy"
        tag { id }

        publishDir "${params.outdir}/variants", mode: "symlink", pattern: "$id"

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

        publishDir "${params.outdir}/variants", mode: "symlink", pattern: "$id"

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

        publishDir "${params.outdir}/alignment", mode: "copy", pattern: "snp.core.vcf"
        publishDir "${params.outdir}/alignment", mode: "copy", pattern: "snp.core.fasta"
        publishDir "${params.outdir}/alignment", mode: "symlink", pattern: "wgs.core.fasta"

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