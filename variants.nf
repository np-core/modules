process EvaluateRandomForest {

    label "forest_evaluate"
    tag { "$id" }

    publishDir "${params.outdir}/forest/evaluation", mode: "copy", pattern: "*.tsv"

    input:
    tuple val(id), file("snippy/*"), file("ont/*"), file("ont/*")
    each model

    output:
    tuple val(id), file("${id}_${model}_application_truth.tsv"), file("${id}_${model}_classifier_truth.tsv"), file("${id}_${model}_${params.eval_caller}_truth.tsv")

    """
    np variants forest-evaluate --prefix ${id}_${model} --snippy_dir snippy/ --ont_dir ont/ --model $model --outdir ${id}_eval --mask_weak $params.eval_mask_weak --caller $params.eval_caller
    mv ${id}_eval/evaluation/${id}_${model}_application_truth.tsv .
    mv ${id}_eval/evaluation/${id}_${model}_classifier_truth.tsv .
    mv ${id}_eval/evaluation/${id}_${model}_${params.eval_caller}_truth.tsv ${id}_${model}_caller_truth.tsv
    """

}

process ProcessRandomForestEvalutions {


    label "forest_evaluate"
    tag { "$id" }

    publishDir "${params.outdir}/forest", mode: "copy", pattern: "rff_application_evaluation.tsv"
    publishDir "${params.outdir}/forest", mode: "copy", pattern: "rff_classfier_evaluation.tsv"
    publishDir "${params.outdir}/forest", mode: "copy", pattern: "rff_caller_evaluation.tsv"

    input:
    tuple val(id), file(app_truth), file(class_truth), file(base_truth)

    output:
    file("*_evaluation.tsv")

    """
    np utils combine-df --dir . --glob "*_application_truth.tsv" --extract "_application_truth.tsv" --extract_split "_" --extract_head "id,model" --output rff_application_evaluation.tsv
    np utils combine-df --dir . --glob "*_classifier_truth.tsv" --extract "_classifier_truth.tsv" --extract_split "_" --extract_head "id,model" --output rff_classfier_evaluation.tsv
    np utils combine-df --dir . --glob "*_caller_truth.tsv" --extract "_caller_truth.tsv" --extract_split "_" --extract_head "id,model" --output rff_caller_evaluation.tsv
    """


}