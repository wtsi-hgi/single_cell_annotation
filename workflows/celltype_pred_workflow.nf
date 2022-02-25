
include { run_keras_celltype_prediction } from '../modules/celltype_prediction/run_keras_cell_type_prediction.nf'

workflow celltype_pred_workflow {
    take:
	keras_input
    main:
        // Identify multiplets using scrublet.
    run_keras_celltype_prediction(keras_input)
    
    emit:
    to_merge = run_keras_celltype_prediction.out.to_merge
        //val(biopsy_type),
        //val(experiment_id),
        //path("*predictions.h5ad"), 
        //emit: to_merge)
}
