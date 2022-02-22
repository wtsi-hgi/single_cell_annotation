
include { run_keras_cell_type_prediction } from '../modules/celltype_prediction/run_keras_cell_type_prediction.nf'

workflow celltype_pred_workflow {
    take:
	keras_input
    main:
        // Identify multiplets using scrublet.
    run_keras_cell_type_prediction(keras_input)
    
    //emit:
    //    multiplet_calls = run_scrublet.out.multiplet_calls
}
