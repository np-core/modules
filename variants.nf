process EvaluateRandomForest {

    label "forest_evaluate"
    tag { "$id" }

    memory { params.forest_evaluate_mem * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 5

    publishDir "${params.outdir}/${ref}/evaluation/${eval_set}/evaluations", mode: "copy", pattern: "*.tsv"

    input:
    tuple val(eval_set), val(id), val(ref), file("snippy/*"), file("ont/*"), file("ont/*")
    each file(model)

    output:
    tuple file("result/evaluation/${id}.${model.simpleName}.${eval_set}.${ref}_application_truth.tsv"), file("result/evaluation/${id}.${model.simpleName}.${eval_set}.${ref}_classifier_truth.tsv"), file("result/evaluation/${id}.${model.simpleName}.${eval_set}.${ref}_${params.caller}_truth.tsv")
    

    """
    np variants forest-evaluate --prefix ${id}.${model.simpleName}.${eval_set}.${ref} --dir_snippy snippy/ --dir_ont ont/ --model $model --mask_weak $params.mask_weak --caller $params.caller --outdir result
    """

}

process ProcessEvaluations {


    label "forest_evaluate"
    tag { "$id" }

    memory { params.forest_evaluate_mem * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.outdir}/", mode: "copy", pattern: "model_evaluation.tsv"
    publishDir "${params.outdir}/", mode: "copy", pattern: "${params.caller}_evaluation.tsv"

    input:
    file(collected)

    output:
    file("*_evaluation.tsv")

    """
    np utils combine-df --dir . --glob "*_application_truth.tsv" --extract "_application_truth.tsv" --extract_split "." --extract_head "id,model,eval_set,reference" --output model_evaluation.tsv
    np utils combine-df --dir . --glob "*_classifier_truth.tsv" --extract "_classifier_truth.tsv" --extract_split "." --extract_head "id,model,eval_set,reference" --output classfier_evaluation.tsv
    np utils combine-df --dir . --glob "*_${params.caller}_truth.tsv" --extract "_${params.caller}_truth.tsv" --extract_split "." --extract_head "id,model,eval_set,reference" --clean --output ${params.caller}_evaluation.tsv
    """

}


process RandomForestTraining {

    label "forest_training"
    tag { "$model - $ref - Composite RFC" }

    memory { params.forest_train_mem * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.outdir}/${ref}/polishers/models", mode: "copy", pattern: "${model}_${ref}.composite.sav"
    publishDir "${params.outdir}/${ref}/polishers/models", mode: "copy", pattern: "${model}_${ref}.qual.sav"
    publishDir "${params.outdir}/${ref}/polishers/models", mode: "copy", pattern: "${model}_${ref}_model"

    input:
    tuple val(model), val(ref), file("ont/*"), file("ont/*"), file("snippy/*")

    output:
    tuple val(model), val(ref), file("${model}_${ref}.composite.sav")
    file("${model}_${ref}_model")

    """
    np variants forest-train --dir_snippy snippy/ --dir_ont ont/ --caller ${params.caller} --prefix ${model}_${ref} --test_size ${params.test_size} --outdir model
    mv model/models/${model}_${ref}.composite.sav ${model}_${ref}.composite.sav 
    mv model/models/${model}_${ref}.qual.sav ${model}_${ref}.qual.sav 
    mv model/training ${model}_${ref}_model
    """

}