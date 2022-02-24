process run_keras_celltype_prediction {
    tag { "${experiment_id}" }

    publishDir  path: "${params.outdir}/celltype_prediction/1_keras/${biopsy_type}/${experiment_id}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(
            val(experiment_id),
            val(biopsy_type),
            path(keras_input_h5ad) // from cellbender into scrublet
        )

    output:
    tuple(
        val(experiment_id),
        val(biopsy_type),
        path("*.png"),
        path("*_gtr_*.h5ad"),// do not capture input Pilot_study_of_dissociation_methods_for_human_gut_tissues7980357-scrublet.h5ad
        path("*predictions.h5ad"), 
        emit: out)
    
    script:
        """
python ${projectDir}/../bin/0057-predict_clusters_keras_model-anndata.py \\
	--h5_anndata \"${keras_input_h5ad}\" \\
	--h5_layer \"${params.celltype_prediction.keras.h5_layer}\" \\
	--keras_model \"${params.celltype_prediction.keras.keras_model}\" \\
	--keras_weights_df \"${params.celltype_prediction.keras.keras_weights_df}\" \\
	--keras_model_cluster_labels \"${params.celltype_prediction.keras.keras_model_cluster_labels}\" \\
	--filter_top_cell_probabilities \"${params.celltype_prediction.keras.filter_top_cell_probabilities}\" \\
	--output_file \"${params.celltype_prediction.keras.output_file}\" ${params.celltype_prediction.keras.save_all_probabilities} 
        """
}
