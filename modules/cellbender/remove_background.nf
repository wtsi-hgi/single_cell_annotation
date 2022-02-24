process remove_background {
  // Remove ambient RNA
  // ------------------------------------------------------------------------
  //tag { output_dir }
  //cache false    // cache results from run
  //maxForks 2   // hard to control memory usage. limit to 3 concurrent
// cb_plot_input,out_paths,results_list,experimentid_outdir_cellbenderunfiltered_expectedcells_totaldropletsinclude

  publishDir  path: "${params.outdir}/cellbender/2_remove_background/${biopsy_type}/${experiment_id}",
      saveAs: {filename ->
        if (filename == "barcodes.tsv.gz") {
        null
        } else if(filename.equalsIgnoreCase("features.tsv.gz")) {
        null
        } else if(filename.equalsIgnoreCase("matrix.mtx.gz")) {
        null
        } else {
        filename.replaceAll("", "")
        }
      },
      mode: "${params.copy_mode}",
      overwrite: "true"

  input:
    tuple(
	val(experiment_id),
	path(file_10x_barcodes),
	path(file_10x_features),
	path(file_10x_matrix),
	path(expected_cells),
	path(total_droplets_include),
        val(biopsy_type),
	val(epochs),
	val(learning_rate),
	val(zdims),
	val(zlayers),
	val(low_count_threshold),
	val(fpr))

  output:
    tuple(
	val(experiment_id),
	val(fpr),
	path("cellbender.h5"),
	path("cellbender_filtered.h5"),
	val(cb_params),
	path(file_10x_barcodes),
	path(file_10x_features),
	path(file_10x_matrix),
	path(expected_cells),
	path(total_droplets_include),
	val(biopsy_type),
	emit: experimentid_cellbender_to_preprocess
    )

    path("${outfile}_cell_barcodes.csv", emit: barcodes)
    
    path("*.h5", emit: filtered_h5s)
    path(file_10x_barcodes, emit: raw_tenx_barcodes)
    path(file_10x_features, emit: raw_tenx_features)
    path(file_10x_matrix, emit: raw_tenx_matrix)
    path("plots/*.png") optional true
    path("plots/*.pdf") optional true


  script:

    lr_string = "${learning_rate}".replaceAll("\\.", "pt")
    lr_string = "${lr_string}".replaceAll("-", "neg")
    fpr_string = "${fpr}".replaceAll("\\.", "pt").replaceAll(" ", "_")
    cb_params = "cellbender_params"
    cb_params = "${cb_params}-epochs_${epochs}"
    cb_params = "${cb_params}__learnrt_${lr_string}"
    cb_params = "${cb_params}__zdim_${zdims}"
    cb_params = "${cb_params}__zlayer_${zlayers}"
    cb_params = "${cb_params}__lowcount_${low_count_threshold}"
    outfile = "cellbender"
    
"""
# LD_PRELOAD to fix mkl/anaconda python error
# cf. https://stackoverflow.com/questions/36659453/intel-mkl-fatal-error-cannot-load-libmkl-avx2-so-or-libmkl-def-so
export LD_PRELOAD=/opt/conda/envs/conda_cellbender/lib/libmkl_core.so:/opt/conda/envs/conda_cellbender/lib/libmkl_sequential.so

rm -fr plots
mkdir txd_input
ln --physical ${file_10x_barcodes} txd_input/barcodes.tsv.gz
ln --physical ${file_10x_features} txd_input/features.tsv.gz
ln --physical ${file_10x_matrix} txd_input/matrix.mtx.gz

cellbender remove-background --input txd_input \\
  --cuda --output ${outfile} \\
  --expected-cells \$(cat ${expected_cells}) \\
  --total-droplets-included \$(cat ${total_droplets_include}) \\
  --model full \\
  --z-dim ${zdims} \\
  --z-layers ${zlayers} \\
  --low-count-threshold ${low_count_threshold} \\
  --epochs ${epochs} \\
  --learning-rate ${learning_rate} \\
  --fpr ${fpr}

# If outfile does not have h5 appended to it, move it.
[ -f ${outfile} ] && mv ${outfile} ${outfile}.h5

mkdir plots
mv *pdf plots/ 2>/dev/null || true
mv *png plots/ 2>/dev/null || true
"""
}
