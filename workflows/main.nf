nextflow.enable.dsl=2

// All input parameters are read from Nextflow config file "inputs.nf"

include { cellbender_workflow as CB_1 } from './cellbender_workflow.nf'
include { multiplet_workflow as MU_2 } from './multiplet_workflow.nf'
include { celltype_pred_workflow as CT_3 } from './celltype_pred_workflow.nf'

workflow {
    
    channel_samples_meta = Channel
	.fromPath(params.samples_metainfo_tsv)
	.splitCsv(header: true, sep: '\t')
	.map{row->tuple(row.sanger_sample_id, row.biopsy_type, row.disease_status)}
	.unique()
    
    channel_samples_meta
	.count().view { "\n --- processing $it samples with biopsy_type and disease_status" }
    
    channel_samples_meta
        .map { a,b,c -> "${a},${b},${c}" }
	.collectFile(name: 'sample_biopsy_disease.csv', newLine: true, sort: true,
		     seed: "sample,biopsy_type,disease_status", storeDir:params.outdir)
        .subscribe { println "\n --- wrote samples with biopsy_type,disease_status to csv table ${it.name}" }

    channel_samples_meta
	.map { sample,biopsy_type,disease_status -> tuple(sample, biopsy_type) }
	.filter { it[1] == "TI" | it[1] == "Rectum" }
	.filter { it[0] != "Not_sequenced" }
        .set {sample_tirectum}
    
    sample_tirectum
        .map { a,b -> "${a},${b}" }
	.collectFile(name: 'sample_biopsy_only_ti_or_rectum.csv', newLine: true, sort: true,
		     seed: "sample,biopsy_type", storeDir:params.outdir)
        .subscribe { println "\n --- wrote Ti and Rectum only samples with biopsy_type,disease_status to csv table ${it.name}" }
    
    sample_tirectum
	.count().view { "\n --- $it samples have TI or Rectum status" }

    sample_tirectum
	.filter { it[1] == "TI" } 
	.map { sample,biopsy_type ->
	tuple(sample, //0
	      biopsy_type, //1
	      params.cellbender.estimate_params_umis.ti.expected_nemptydroplets_umi_cutoff, //2
	      params.cellbender.estimate_params_umis.ti.method_estimate_ncells,//3
	      params.cellbender.estimate_params_umis.ti.lower_bound_umis_estimate_ncells,//4
	      params.cellbender.estimate_params_umis.ti.method_estimate_nemptydroplets,//5
	      params.cellbender.estimate_params_umis.ti.lower_bound_umis_estimate_nemptydroplets,//6
	      params.cellbender.estimate_params_umis.ti.upper_bound_umis_estimate_nemptydroplets,//7
	      params.cellbender.estimate_params_umis.ti.estimate_nemptydroplets_umi_add_factor,//8
	      params.cellbender.estimate_params_umis.ti.estimate_nemptydroplets_subtract_cell_factor,//9
	      params.cellbender.estimate_params_umis.ti.estimate_nemptydroplets_min_drop,//10
	      params.cellbender.estimate_params_umis.ti.epochs,//11
	      params.cellbender.estimate_params_umis.ti.learning_rate,//12
	      params.cellbender.estimate_params_umis.ti.zdim,//13
	      params.cellbender.estimate_params_umis.ti.zlayers,//14
	      params.cellbender.estimate_params_umis.ti.low_count_threshold,//15
	      params.cellbender.estimate_params_umis.ti.fpr) }//16
	.set {sample_ti_cellbender_params}

    sample_tirectum
	.filter { it[1] == "TI" } 
	.map { sample,biopsy_type ->
	tuple(sample, biopsy_type,
	      params.multiplet.scrublet.ti.expected_multiplet_rate,//2
	      params.multiplet.scrublet.ti.n_simulated_multiplet,//3
	      params.multiplet.scrublet.ti.multiplet_threshold_method,//4
	      params.multiplet.scrublet.ti.scale_log10//5
	      ) }
	.set {sample_ti_scrublet_params}

    sample_tirectum
	.filter { it[1] == "Rectum" } 
	.map { sample,biopsy_type ->
	tuple(sample, biopsy_type,
	      params.cellbender.estimate_params_umis.rectum.expected_nemptydroplets_umi_cutoff,
	      params.cellbender.estimate_params_umis.rectum.method_estimate_ncells,
	      params.cellbender.estimate_params_umis.rectum.lower_bound_umis_estimate_ncells,
	      params.cellbender.estimate_params_umis.rectum.method_estimate_nemptydroplets,
	      params.cellbender.estimate_params_umis.rectum.lower_bound_umis_estimate_nemptydroplets,
	      params.cellbender.estimate_params_umis.rectum.upper_bound_umis_estimate_nemptydroplets,
	      params.cellbender.estimate_params_umis.rectum.estimate_nemptydroplets_umi_add_factor,
	      params.cellbender.estimate_params_umis.rectum.estimate_nemptydroplets_subtract_cell_factor,
	      params.cellbender.estimate_params_umis.rectum.estimate_nemptydroplets_min_drop,
	      params.cellbender.estimate_params_umis.rectum.epochs,
	      params.cellbender.estimate_params_umis.rectum.learning_rate,
	      params.cellbender.estimate_params_umis.rectum.zdim,
	      params.cellbender.estimate_params_umis.rectum.zlayers,
	      params.cellbender.estimate_params_umis.rectum.low_count_threshold,
	      params.cellbender.estimate_params_umis.rectum.fpr) }
	.set {sample_rectum_cellbender_params}

    sample_tirectum
	.filter { it[1] == "Rectum" } 
	.map { sample,biopsy_type ->
	tuple(sample, biopsy_type,
	      params.multiplet.scrublet.rectum.expected_multiplet_rate,
	      params.multiplet.scrublet.rectum.n_simulated_multiplet,
	      params.multiplet.scrublet.rectum.multiplet_threshold_method,
	      params.multiplet.scrublet.rectum.scale_log10
	      ) }
	.set {sample_rectum_scrublet_params}

    sample_ti_cellbender_params.mix(
	sample_rectum_cellbender_params).set {
	sample_tirectum_cellbender_params}

    // sample_tirectum_cellbender_params.view()

    sample_ti_scrublet_params.mix(
	sample_rectum_scrublet_params).set {
	sample_tirectum_scrublet_params}

    // sample_tirectum_scrublet_params.view()
    
    channel_cellranger_dirs = Channel
	.fromPath("${params.cellranger_10x_dir}/**", maxDepth: 3)
	.map { dir_ -> file(dir_).getParent() } // list cellranger dirs pulled from Irods
	.map { dir_ -> [file(dir_).getName(), file(dir_)] } // attach samplename
	.filter {it[0] != "full_data" }
	.map { sample, dir_ -> [sample,
                                // attach filtered and raw cellranger subdirs
                                //  checkIfExists: true to exit if file not found!
				file("${dir_}/raw_feature_bc_matrix", checkIfExists: true),
				file("${dir_}/filtered_feature_bc_matrix", checkIfExists: true),
				file("${dir_}/metrics_summary.csv", checkIfExists: true) ] }
	.unique()
        .join(sample_tirectum)
	.map {sample, raw, filtered, metrics, biopsy_type -> tuple(sample,raw,filtered,metrics)}

    //channel_cellranger_dirs
    //    .view()
    
    // extract cellranger estimated number of cells
    // which is file column of metrics_summary.csv file:
    channel_cellranger_dirs
        .map { a,b,c,d -> [a, file("${d}").readLines()] } // read lines of metrics_summary.csv
        .map { a,b -> [a, b[1]] } // keep only second row (first one is header)
        .map { row -> [row[0],row[1].replaceAll(/",.*$/,'')] } // keep only first csv column, which is Estimated Number of Cells
        .map { row -> [row[0],row[1].replaceAll(/^"/,'')] } // remove extra "
        .map { row -> [row[0],row[1].replaceAll(/,/,'')] } // remove extra ,
        .map { row -> [row[0],row[1].replaceAll(/".*/,'')] } // remove extra 2nd column when n cells < 1000
	.take(4).set {sample_cellranger_estimated_n_cells}
    
    sample_cellranger_estimated_n_cells
	.count().view { "\n --- extracted $it n cellranger_estimated_n_cells from metrics_summary.csv files " }
    
    sample_cellranger_estimated_n_cells
        .map { a,b-> "${a},${b}" }
	.collectFile(name: 'sample_cellranger_estimated_n_cells.csv', newLine: true, sort: true,
		     seed: "sample,cellranger_estimated_n_cells", storeDir:params.outdir)
        .subscribe { println "\n --- wrote samples cellranger_estimated_n_cells to csv table ${it.name}" }

    channel_cellranger_dirs
	.count().view { "\n --- processing $it samples cellranger directories" }
    
    channel_cellranger_dirs
        .map { a,b,c,d -> "${a},${b},${c},${d}" }
	.collectFile(name: 'sample_raw_feature_metrics.csv', newLine: true, sort: true,
		     seed: "sample,raw_feature_bc_matrix,filtered_feature_bc_matrix,metrics_summary", storeDir:params.outdir)
        .subscribe { println "\n --- wrote samples cellranger directories to csv table ${it.name}" }

	channel_cellranger_dirs
            .map { sample ,raw ,c ,d -> [sample,
					 file("${raw}/barcodes.tsv.gz", checkIfExists: true ),
					 file("${raw}/features.tsv.gz", checkIfExists: true),
					 file("${raw}/matrix.mtx.gz", checkIfExists: true)] }
	    .set{ channel__file_paths_10x }

	channel_cellranger_dirs
            .map { sample ,raw ,c ,d -> [sample,raw] }
	    .set{ ch_experimentid_paths10x_raw }

	channel_cellranger_dirs
            .map { sample ,raw , filt ,d -> [sample,filt] }
	    .set{ ch_experimentid_paths10x_filtered }

	// here pass in the number of cells detected by cellranger/ 
    if (params.cellbender.estimate_params_umis.use_cellranger_estimate){
	    log.info '\n --- use_cellranger_estimate=true -> use cellranger estimated number of cells'
            sample_cellranger_estimated_n_cells.set{ ncells_cellranger }}
	else {
	    log.info '\n --- use_cellranger_estimate=false -> not using cellranger estimated number of cells'
	    sample_cellranger_estimated_n_cells.map{ a,b -> [a,'0']}.set{ ncells_cellranger } }


    // all prepared channels:
//    ch_experimentid_paths10x_raw.take(1).view()
//    ch_experimentid_paths10x_filtered.take(1).view()
//    channel__file_paths_10x.take(1).view()
//    ncells_cellranger.take(1).view()
//    sample_tirectum_cellbender_params.take(1).view()
//    sample_tirectum_scrublet_params.take(1).view()
    
    if (params.run_cellbender_workflow) {
	log.info '\n ---- running cellbender workflow ---- '
	CB_1(
	    // sample, path to cellranger raw_feature_bc_matrix folder
	    ch_experimentid_paths10x_raw,
	    
            //sample, filtered_feature_bc_matrix 
            ch_experimentid_paths10x_filtered,
	    
	    // sample, path to raw_feature_bc_matrix/barcodes.tsv.gz
	    //       , raw_feature_bc_matrix/features.tsv.gz
	    //       , raw_feature_bc_matrix/matrix.mtx.gz
	    channel__file_paths_10x,
	    
	    // row.experiment_id, row.Estimated_Number_of_Cells
	    // sample, cellranger n cells or '0'
	    ncells_cellranger,
	    
	    // sample, cellbender param 1 ... n
	    sample_tirectum_cellbender_params) 

	if (params.run_multiplet_workflow) {
	    log.info '\n --- running multiplet workflow ---- '
            MU_2(CB_1.out.filt10x
		    .join(sample_tirectum_scrublet_params, by: [0,1]))

	    if (params.run_celltype_pred_workflow) {
		log.info '\n ---- running keras cell type prediction workflow ---- '
	  	CT_3(MU_2.out.keras_input)}}}}

workflow.onError {
    log.info "\n\n\n --- Pipeline execution stopped with the following message: ${workflow.errorMessage}" }

workflow.onComplete {
    log.info "\n --- Pipeline completed at: $workflow.complete"
    log.info "\n --- Command line: $workflow.commandLine"
    log.info "\n --- Execution status: ${ workflow.success ? 'OK' : 'failed' }\n"

    if (params.on_complete.sync_results_to_gitlab && ("${workflow.success}" == "OK")) {
	log.info "\n --- Sync results to gitlab: run script bin/sync_results_to_gitlab.sh\n"
	def proc = "bash ${projectDir}/../bin/sync_results_to_gitlab.sh ${projectDir}/../results".execute()
	def b = new StringBuffer()
	proc.consumeProcessErrorStream(b)
	log.info proc.text
	log.info b.toString() }}
