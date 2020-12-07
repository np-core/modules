process EvaluateRandomForest {

    label "forest_evaluate"
    tag { "$id" }

    memory { params.forest_evaluate_mem * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 5

    publishDir "${params.outdir}/evaluations/data", mode: "copy", pattern: "*.tsv"

    input:
    tuple val(id), file("snippy/*"), file("ont/*"), file("ont/*")
    each file(model)

    output:
    tuple file("${id}.${model.simpleName}.application.truth.tsv"), file("${id}.${model.simpleName}.classifier.truth.tsv"), file("${id}.${model.simpleName}.raw.truth.tsv")

    """
    np variants forest-evaluate --prefix ${id}_${model.simpleName} --dir_snippy snippy/ --dir_ont ont/ --model $model --mask_weak $params.mask_weak --caller $params.caller
    """

}

process ProcessRandomForestEvaluations {


    label "forest_evaluate"
    tag { "$id" }

    memory { params.forest_evaluate_mem * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.outdir}/evaluations", mode: "copy", pattern: "rff_application_evaluation.tsv"
    publishDir "${params.outdir}/evaluations/classifier", mode: "copy", pattern: "rff_classfier_evaluation.tsv"
    publishDir "${params.outdir}/evaluations/raw", mode: "copy", pattern: "rff_${params.eval_caller}_evaluation.tsv"

    input:
    file(collected)

    output:
    file("*_evaluation.tsv")

    """
    np utils combine-df --dir . --glob "*.application.truth.tsv" --extract ".application.truth.tsv" --extract_split "." --extract_head "id,model" --output rff_application_evaluation.tsv
    np utils combine-df --dir . --glob "*.classifier.truth.tsv" --extract ".classifier.truth.tsv" --extract_split "." --extract_head "id,model" --output rff_classfier_evaluation.tsv
    np utils combine-df --dir . --glob "*.caller.truth.tsv" --extract ".caller.truth.tsv" --extract_split "." --extract_head "id,model" --clean --output rff_${params.eval_caller}_evaluation.tsv
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