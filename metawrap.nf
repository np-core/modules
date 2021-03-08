process MetaWrap::QC {

    label "metawrap_qc"
    tag { "$id" }

    publishDir "${params.outdir}/metawrap", mode: "symlink", pattern: "READ_QC"

    input:
    tuple val(id), file(fwd), file(rev)

    output:
    tuple val(id), file("READ_QC/final_pure_reads_1.fastq"), file("READ_QC/final_pure_reads_2.fastq")

    """
    metawrap read_qc -1 $fwd -2 $rev -t $task.cpus -o READ_QC/ $params.qc_options
    """

}

process MetaWrap::ASSEMBLY {

    label "metawrap_assembly"
    tag { "$id" }

    publishDir "${params.outdir}/metawrap", mode: "symlink", pattern: "ASSEMBLY"

    input:
    tuple val(id), file(fwd), file(rev)

    output:
    tuple val(id), file("ASSEMBLY/final_assembly.fasta")

    """
    metawrap assembly -1 $fwd -2 $rev-m $task.memory -t $task.cpus $params.assembly_options -o ASSEMBLY
    """

}

process MetaWrap::BINNING {

    label "metawrap_binning"
    tag { "$id" }

    publishDir "${params.outdir}/metawrap", mode: "symlink", pattern: "INITIAL_BINNING"
    publishDir "${params.outdir}/metawrap", mode: "symlink", pattern: "BIN_REFINEMENT"

    input:
    tuple val(id), file(assembly)
    tuple val(id), file(fwd), file(rev)

    output:
    tuple val(id), file("BIN_REFINEMENT")

    """
    metawrap binning -o INITIAL_BINNING -t $task.cpus -a $assembly --metabat2 --maxbin2 --concoct $fwd $rev
    metawrap bin_refinement -o BIN_REFINEMENT -t $task.cpus -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins/ -C INITIAL_BINNING/concoct_bins/ -c $params.completeness -x $params.contamination
    """

}

process MetaWrap::BINASSEMBLY {

    label "metawrap_binassembly"
    tag { "$id" }

    publishDir "${params.outdir}/metawrap", mode: "symlink", pattern: "BIN_REASSEMBLY"


    input:
    tuple val(id), file(bin_refinement)
    tuple val(id), file(fwd), file(rev)

    output:
    tuple val(id), file("BIN_REASSEMBLY/")

    """
    metawrap reassemble_bins -o BIN_REASSEMBLY -1 $fwd -2 $rev -t $task.cpus -m $task.memory -c $params.completeness -x $params.contamination -b $bin_refinement/metawrap_${params.completeness}_${params.contamination}_bins
    """

}

process MetaWrap::BINOPS {

    label "metawrap_binops"
    tag { "$id" }

    publishDir "${params.outdir}/metawrap", mode: "symlink", pattern: "FUNCT_ANNOT"
    publishDir "${params.outdir}/metawrap", mode: "symlink", pattern: "INITIAL_BINNING"


    input:
    tuple val(id), file(bin_reassembly)
    tuple val(id), file(fwd), file(rev)

    output:
    tuple val(id), file

    """
    metawrap classify_bins -b $bin_reassembly/reassembled_bins -o BIN_CLASSIFICATION -t $task.cpus
    metawrap annotate_bins -o FUNCT_ANNOT -t $task.cpus -b $bin_reassembly/reassembled_bins/
    """

}