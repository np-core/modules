process EvaluateRandomForest {

    label "forest_evaluate"
    tag { "$id" }

    publishDir "${params.outdir}/forest/evaluation", mode: "copy", pattern: "*.tsv"

    input:
    tuple val(id), file("snippy/*"), file("ont/*"), file("ont/*")
    each file(model)

    output:
    tuple file("${id}.${model.simpleName}.application.truth.tsv"), file("${id}.${model.simpleName}.classifier.truth.tsv"), file("${id}.${model.simpleName}.caller.truth.tsv")

    """
    np variants forest-evaluate --prefix ${id}_${model.simpleName} --dir_snippy snippy/ --dir_ont ont/ --model $model --outdir ${id}_eval --mask_weak $params.eval_mask_weak --caller $params.eval_caller
    mv ${id}_eval/evaluation/${id}_${model.simpleName}_application_truth.tsv ${id}.${model.simpleName}.application.truth.tsv
    mv ${id}_eval/evaluation/${id}_${model.simpleName}_classifier_truth.tsv ${id}.${model.simpleName}.classifier.truth.tsv
    mv ${id}_eval/evaluation/${id}_${model.simpleName}_${params.eval_caller}_truth.tsv ${id}.${model.simpleName}.caller.truth.tsv
    """

}

process ProcessRandomForestEvaluations {


    label "forest_evaluate"
    tag { "$id" }

    publishDir "${params.outdir}/forest", mode: "copy", pattern: "rff_application_evaluation.tsv"
    publishDir "${params.outdir}/forest", mode: "copy", pattern: "rff_classfier_evaluation.tsv"
    publishDir "${params.outdir}/forest", mode: "copy", pattern: "rff_${params.eval_caller}_evaluation.tsv"

    input:
    file(collected)

    output:
    file("*_evaluation.tsv")

    """
    np utils combine-df --dir . --glob "*.application.truth.tsv" --extract ".application.truth.tsv" --extract_split "." --extract_head "id,model" --output rff_application_evaluation.tsv
    np utils combine-df --dir . --glob "*.classifier.truth.tsv" --extract ".classifier.truth.tsv" --extract_split "." --extract_head "id,model" --output rff_classfier_evaluation.tsv
    np utils combine-df --dir . --glob "*.caller.truth.tsv" --extract ".caller.truth.tsv" --extract_split "." --extract_head "id,model" --output rff_${params.eval_caller}_evaluation.tsv
    """

}


process TrainRandomForest {

    label "forest_training"
    tag { "$model_name" }

    publishDir "${params.outdir}/forest/training", mode: "copy", pattern: "${model_name}_${reference_name}_model"

    input:
    tuple val(model_name), val(reference_name), file("ont/*"), file("ont/*"), file("snippy/*")

    output:
    file("*.sav")

    """
    np variants forest-train --dir_snippy snippy/ --dir_ont ont/ --caller ${params.caller} --prefix ${model_name}_${reference_name} --test_size ${params.test_size} --outdir ${model_name}_${reference_name}_model
    """

}