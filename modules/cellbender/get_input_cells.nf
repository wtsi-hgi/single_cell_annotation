process get_input_cells {
  // Calculates thresholds for input cells of remove_background
  // ------------------------------------------------------------------------

    publishDir path: "${params.outdir}/cellbender/1_get_input_cells/${biopsy_type}/${experiment_id}",
	mode: "${params.copy_mode}",
	overwrite: "true"
    
    input:
    tuple(val(experiment_id),
          path(file_10x_barcodes),
          path(file_10x_features),
          path(file_10x_matrix),
	  val(ncells),
          val(biopsy_type),
	  val(method_estimate_ncells),
	  val(expected_nemptydroplets_umi_cutoff),
	  val(method_estimate_ncells),
	  val(lower_bound_umis_estimate_ncells),
	  val(method_estimate_nemptydroplets),
	  val(lower_bound_umis_estimate_nemptydroplets),
	  val(upper_bound_umis_estimate_nemptydroplets),
	  val(estimate_nemptydroplets_umi_add_factor),
	  val(estimate_nemptydroplets_subtract_cell_factor),
	  val(estimate_nemptydroplets_min_drop))

  output:
    tuple(
    val(experiment_id),
    path(file_10x_barcodes),
    path(file_10x_features),
    path(file_10x_matrix),
    path("${outfile}-expected_cells.txt"),
    path("${outfile}-total_droplets_included.txt"),
    emit: cb_input
    )
    path(
    "${outfile}-expected_cells.txt",
    emit: expected_cells
    )
    path(
    "${outfile}-total_droplets_included.txt",
    emit: total_droplets_include
    )
    path("${outfile}-cell_estimate_cutoff.tsv.gz")
    path("${outfile}-total_droplets_cutoff.tsv.gz")
    path("plots/*.png") optional true
    path("plots/*.pdf") optional true

  script:
    outfile = "umi_count_estimates"
    cmd__expected_ncells = ""
    cmd__droplets_include = ""

    if (method_estimate_ncells=='expected'){
        ncells2=ncells.toInteger()-1000
        cell_numbers = "--expected_ncells ${ncells2}"
    }else{
        cell_numbers = ""
    }

    """

    rm -fr plots
    mkdir txd_input
    ln --physical ${file_10x_barcodes} txd_input/barcodes.tsv.gz
    ln --physical ${file_10x_features} txd_input/features.tsv.gz
    ln --physical ${file_10x_matrix} txd_input/matrix.mtx.gz

   python ${projectDir}/../bin/get_estimates_from_umi_counts.py ${cell_numbers} \\
     --tenxdata_path txd_input \\
     --output_file ${outfile} \\
     --expected_nemptydroplets_umi_cutoff ${expected_nemptydroplets_umi_cutoff} \\
     --method_estimate_ncells ${method_estimate_ncells} \\
     --lower_bound_umis_estimate_ncells ${lower_bound_umis_estimate_ncells} \\
     --method_estimate_nemptydroplets ${method_estimate_nemptydroplets} \\
     --lower_bound_umis_estimate_nemptydroplets ${lower_bound_umis_estimate_nemptydroplets} \\
     --upper_bound_umis_estimate_nemptydroplets ${upper_bound_umis_estimate_nemptydroplets} \\
     --estimate_nemptydroplets_add_umifactor ${estimate_nemptydroplets_umi_add_factor} \\
     --estimate_nemptydroplets_subtract_dropletfactor ${estimate_nemptydroplets_subtract_cell_factor} \\
     --estimate_nemptydroplets_min_nemptydroplets ${estimate_nemptydroplets_min_drop} ${cmd__expected_ncells} ${cmd__droplets_include}

    mkdir plots
    mv *pdf plots/ 2>/dev/null || true
    mv *png plots/ 2>/dev/null || true
    """
}

