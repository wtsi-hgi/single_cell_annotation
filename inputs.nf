params {

    // path to directory where cellranger 10x directories are 3 depths down, where name of 10x dir is name of sample.
    cellranger_10x_dir = '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/results/iget_cellranger/full_data'

    // path to tab-delimited table with columns row.sanger_sample_id, row.biopsy_type, row.disease_status
    // For now, processes samples where row.biopsy_type is "Ti" or "Rectum". Ignores other samples
    // row.sanger_sample_id is sample name which must match a 10x dir name.
    samples_metainfo_tsv = '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/results/info_stats/samples_metainfo.tsv'

    // which subworkflows to run:
    run_cellbender_workflow = true
    run_multiplet_workflow = true // can only run if run_cellbender_workflow has run
    run_celltype_pred_workflow = true // can only run if run_multiplet_workflow has run
    on_complete { sync_results_to_gitlab = true } // push some outputs to gitlab repo in workflow.onComplete{}
    
    copy_mode = "rellink" // publish relative links from work to results publish dir
    mem_gpu = 10000 // Gb RAM memory for cellbender gpu task

    cellbender {
	estimate_params_umis {
	    use_cellranger_estimate = false // TODO: set to true to use cellranger n estimated cells

	    // "Ti" or "Rectum" biopsy_type samples are processed with different parameters:
	    rectum {
                // input descriptions at
		// https://github.com/andersonlab/letaylor/tree/master/studies/single_cell_continuous_data_integration
		expected_nemptydroplets_umi_cutoff = 0
		method_estimate_ncells = 'dropletutils::barcoderanks::inflection'
		//method_estimate_ncells = 'cellrangerv2::expected' //this method feeds in the cellranger estimate ncells
		//method_estimate_ncells ='expected'
		lower_bound_umis_estimate_ncells = 1000
		method_estimate_nemptydroplets = 'dropletutils::barcoderanks::inflection,dropletutils::barcoderanks::knee,0.33'
		lower_bound_umis_estimate_nemptydroplets = 10
		upper_bound_umis_estimate_nemptydroplets = 250
		estimate_nemptydroplets_umi_add_factor = 0
		estimate_nemptydroplets_subtract_cell_factor = 0
		estimate_nemptydroplets_min_drop = 0
		epochs = 250
		learning_rate = 0.0000001
		zdim = 100
		zlayers = 500
		low_count_threshold = 10
		fpr = "0.1"} // value="0.01 0.05 0.1"

	    ti {
		expected_nemptydroplets_umi_cutoff = 0
		method_estimate_ncells = 'dropletutils::barcoderanks::inflection'
		//method_estimate_ncells = 'cellrangerv2::expected' //this method feeds in the cellranger estimate ncells
		//method_estimate_ncells ='expected'
		lower_bound_umis_estimate_ncells = 1000
		method_estimate_nemptydroplets = 'dropletutils::barcoderanks::inflection,dropletutils::barcoderanks::knee,0.33'
		lower_bound_umis_estimate_nemptydroplets = 10
		upper_bound_umis_estimate_nemptydroplets = 250
		estimate_nemptydroplets_umi_add_factor = 0
		estimate_nemptydroplets_subtract_cell_factor = 0
		estimate_nemptydroplets_min_drop = 0
		epochs = 250
		learning_rate = 0.0000001
		zdim = 100
		zlayers = 500
		low_count_threshold = 10
		fpr = "0.1"} // value="0.01 0.05 0.1"
	}
    }
    multiplet {
	scrublet {
	    // "Ti" or "Rectum" biopsy_type samples are processed with different parameters:
	    rectum {
		expected_multiplet_rate = 0.1
		n_simulated_multiplet = 100000
		multiplet_threshold_method = 'threshold_li'
		scale_log10 = "False" }
	    ti {
		expected_multiplet_rate = 0.1
		n_simulated_multiplet = 100000
		multiplet_threshold_method = 'threshold_li'
		scale_log10 = "False" }}
	
    }

    celltype_prediction {

	// parameters for the keras celltype prediction python script
        keras {
	    keras_model = '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy/sc_qc_cluster/pca_plus5/cellbender_fpr0pt1_filteroutlier_0pt1/nf_results/normalize=total_count.vars_to_regress=none.hvg_exclude=data-variable_gene_filter.scores=data-gene_scores/reduced_dims-vars_to_regress=none-bbknn.batch=experiment_id.n_pcs=29/cluster.number_neighbors=-1.method=leiden.resolution=3.0/validate_resolution/adata-normalized_pca-bbknn-umap-clustered-sparsity_l1=1pt0E-4-train_size_cells=-1.h5'
	    keras_weights_df = '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy/sc_qc_cluster/pca_plus5/cellbender_fpr0pt1_filteroutlier_0pt1/nf_results/normalize=total_count.vars_to_regress=none.hvg_exclude=data-variable_gene_filter.scores=data-gene_scores/reduced_dims-vars_to_regress=none-bbknn.batch=experiment_id.n_pcs=29/cluster.number_neighbors=-1.method=leiden.resolution=3.0/validate_resolution/adata-normalized_pca-bbknn-umap-clustered-sparsity_l1=1pt0E-4-train_size_cells=-1-weights.tsv.gz'
            h5_layer = 'log1p_cp10k'
	    // cloned from github:
	    // https://github.com/andersonlab/sc_nextflow-studies/blob/main/gut-freeze003/ti-cd_healthy/results/cluster_annotations/data-cluster_labels.csv
	    keras_model_cluster_labels = '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/single_cell_annotation/inputs/sc_nextflow-studies-main/gut-freeze003/ti-cd_healthy/results/cluster_annotations/data-cluster_labels.csv'
	    filter_top_cell_probabilities = '0.5,0.75'
	    output_file = 'cellbender_fpr0pt1-scrublet-ti_freeze003_prediction'
	    save_all_probabilities = '--save_all_probabilities'
	}
    }

}
