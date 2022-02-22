process run_keras_cell_type_prediction {
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
        val(experiment_id, emit: experiment_id)
        path("${experiment_id}-scrublet.tsv.gz", emit: multiplet_calls)
        path(
            "${experiment_id}-multiplet_calls_published.txt",
            emit: multiplet_calls_published
        )

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
