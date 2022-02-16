params {
    cellranger_10x_dir = '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/results/iget_cellranger/full_data'
    samples_metainfo_tsv = '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/results/info_stats/samples_metainfo.tsv'
    copy_mode = "rellink"
    mem_gpu = 10000
    run_cellbender_workflow = true
    run_multiplet_workflow = true // can only run if run_cellbender_workflow has run
    run_celltype_pred_workflow = false // can only run if run_multiplet_workflow has run

    cellbender {
	description= 'Parameters for cellbender remove background.'
	estimate_params_umis {
	    
	    method_estimate_ncells = '0' // set to 'expected' to use cellranger n estimated cells
            description = """
                Cellbender requires two parameters=
                    (1) expected-cells
                    (2) total-droplets-included
                Expected number of cells: expected a priori from the experimental
                    design or based on the UMI curve at a point where one is
                    reasonably sure that all droplets to the left on the UMI
                    curve are real cells.
                Total droplets included: emtpy droplets. Point on the UMI curve
                    where every droplet to the right of this number on
                    the UMI curve should be surely empty.
                There are several ways to provide these parameters to this
                workflow:
                    (1) In the file_paths_10x.tsv file under the ncells_expected
                        ndroplets_include_cellbender columns
                    (2) Estimate both parameters via the
                        get_estimates_from_umi_counts.py script.
                get_estimates_from_umi_counts.py options for # of cells:
                    (1) method_estimate_ncells: method used to estimate knee or
                        inflection point on UMI plot.
                    (2) lower_bound_umis_estimate_ncells: remove cells with UMIs
                        below this bound before estimating the knee/inflection.
                get_estimates_from_umi_counts.py options for # of empty droplets:
                    (1) method_estimate_ncells: method used to estimate knee or
                        inflection point on UMI plot.
                    (2) lower_bound_umis_estimate_ncells: remove cells with UMIs
                        below this bound before estimating the knee/inflection.
                    (3) upper_bound_umis_estimate_ncells: remove cells with UMIs
                        above this bound before estimating the knee/inflection.
                    (4) estimate_nemptydroplets_umi_add_factor: after
                        identifying an inflection point add this number to
                        the umi counts to get the final cutoff point,
                    (5) expected_nemptydroplets_umi_cutoff: set the empty droplet
                        value using this umi cutoff.
                    (6) estimate_nemptydroplets_min_drop: if the estimated droplet
                        cutoff is < this value, then set it to this value.
                    (7) estimate_nemptydroplets_subtract_cell_factor: subtract this
                        number from the total number of cell estimates.
                If exact number per sample is provided in file_paths_10x.tsv via
                the ndroplets_include_cellbender column, then estimates are not
                performed and this parameter does nothing.
                NOTE: Setting this parameter will likely vary according
                    tissue / dataset - evaluate the UMI count plots across the
                    full dataset before setting.
                WARNING Do not attempt to iterate over many parameters here as
                    these settings are not recorded in the output dir.
                """

            value{
                expected_nemptydroplets_umi_cutoff= 0
                method_estimate_ncells= 'dropletutils::barcoderanks::inflection'
                //method_estimate_ncells= 'cellrangerv2::expected' //this method feeds in the cellranger estimate ncells
                //method_estimate_ncells='expected'
                lower_bound_umis_estimate_ncells= 1000
                method_estimate_nemptydroplets= 'dropletutils::barcoderanks::inflection,dropletutils::barcoderanks::knee,0.33'
                lower_bound_umis_estimate_nemptydroplets= 10
                upper_bound_umis_estimate_nemptydroplets= 250
                estimate_nemptydroplets_umi_add_factor= 0
                estimate_nemptydroplets_subtract_cell_factor= 0
                estimate_nemptydroplets_min_drop= 0
            }
        }
	
        epochs{
            description = """CellBender parameter. Number of epochs for training.
                    CellBender default is 150."""
            value= 250
        }
        learning_rate{
            description = """CellBender parameter. Learning rate. If lower learning
                    rate, may need to increase the number of epochs.
                    CellBender default is 0.0001."""
            value= 0.0000001
        }
        zdim{
            description= """Dimension of latent variable z, in v2 this parameter
                    influences the prior on cell counts.
                    https://github.com/broadinstitute/CellBender/issues/42
                    CellBender default is 100."""
            value = 100
        }
        zlayers{
            description= """Dimension of hidden layers in the encoder for z.
                    CellBender default is 500."""
            value= 500
        }
        low_count_threshold{
            description="""Droplets with UMI counts below this number are completely
                    excluded from the analysis. This can help identify the correct
                    prior for empty droplet counts in the rare case where empty
                    counts are extremely high (over 200).
                    CellBender default is 15."""
            value=10
        }
        fpr{
            description = """CellBender parameter. A value of 0.01 is generally
                    quite good, but you can generate a few output count matrices
                    and compare them by choosing a few values: 0.01 0.05 0.1.
                    Target false positive rate in (0, 1). A false positive is a true
                    signal count that is erroneously removed. More background removal
                    is accompanied by more signal removal at high values of FPR.
                    You can specify multiple values, which will create multiple
                    output files.
                    CellBender default is 0.01."""
            value="0.1" // value="0.01 0.05 0.1"
            }
	
    }
    
    multiplet {
	description= 'multiplet scrubet parameters'
	scrublet {
	    expected_multiplet_rate = 0.1
	    n_simulated_multiplet = 100000
	    multiplet_threshold_method = 'threshold_li'
	    scale_log10 = "False" // set to string "True" or "False"
	}
    }

    celltype_ {
	description= 'keras'
        keras_model = '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy/sc_qc_cluster/pca_plus5/cellbender_fpr0pt1_filteroutlier_0pt1/nf_results/normalize=total_count.vars_to_regress=none.hvg_exclude=data-variable_gene_filter.scores=data-gene_scores/reduced_dims-vars_to_regress=none-bbknn.batch=experiment_id.n_pcs=29/cluster.number_neighbors=-1.method=leiden.resolution=3.0/validate_resolution/adata-normalized_pca-bbknn-umap-clustered-sparsity_l1=1pt0E-4-train_size_cells=-1.h5'
	keras_weights_df = '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy/sc_qc_cluster/pca_plus5/cellbender_fpr0pt1_filteroutlier_0pt1/nf_results/normalize=total_count.vars_to_regress=none.hvg_exclude=data-variable_gene_filter.scores=data-gene_scores/reduced_dims-vars_to_regress=none-bbknn.batch=experiment_id.n_pcs=29/cluster.number_neighbors=-1.method=leiden.resolution=3.0/validate_resolution/adata-normalized_pca-bbknn-umap-clustered-sparsity_l1=1pt0E-4-train_size_cells=-1-weights.tsv.gz'
	keras_model_cluster_labels = '/nfs/users/nfs_l/lt9/repo/sc_nextflow-studies/gut-freeze003/ti-cd_healthy/results/cluster_annotations/data-cluster_labels.csv'
    }
}
