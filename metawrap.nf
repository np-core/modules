process MetaWrap {

    label "metawrap"
    tag { "$id" }

    publishDir "${params.outdir}/metawrap", mode: "symlink", pattern: "FUNCT_ANNOT"
    publishDir "${params.outdir}/metawrap", mode: "symlink", pattern: "BIN_CLASSIFICATION"
    publishDir "${params.outdir}/metawrap", mode: "symlink", pattern: "BIN_REFINEMENT"
    publishDir "${params.outdir}/metawrap", mode: "symlink", pattern: "INITIAL_BINNING"

    input:
    tuple val(id), file(fwd), file(rev)

    """
    metawrap assembly -1 $fwd -2 $rev -m $param. -t $task.cpus --use-metaspades -o ASSEMBLY
    metawrap binning -o INITIAL_BINNING -t $task.cpus -a ASSEMBLY/final_assembly.fasta --metabat2 --maxbin2 --concoct $fwd $rev
    metawrap bin_refinement -o BIN_REFINEMENT -t $task.cpus -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins/ -C INITIAL_BINNING/concoct_bins/ -c 70 -x 5
    metawrap classify_bins -b BIN_REASSEMBLY/reassembled_bins -o BIN_CLASSIFICATION -t $task.cpus
    metaWRAP annotate_bins -o FUNCT_ANNOT -t $task.cpus -b BIN_REASSEMBLY/reassembled_bins/
    """

}