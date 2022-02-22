
process cellbender__preprocess_output {
    publishDir  path: "${params.outdir}/cellbender/3_preprocess_output/${biopsy_type}/${experiment_id}",
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
	val(fpr),
	path(cellbender_h5),
	path(cellbender_filtered_h5),
	val(cb_params),
	path(file_10x_barcodes),
	path(file_10x_features),
	path(file_10x_matrix),
	path(expected_cells),
	path(total_droplets_include),
	val(biopsy_type))
    
    output:
    tuple(
	val(experiment_id),
	val(biopsy_type),
	path("*filtered_10x_mtx/barcodes.tsv.gz"),
	path("*filtered_10x_mtx/features.tsv.gz"),
	path("*filtered_10x_mtx/matrix.mtx.gz"),
      emit: filt10x
    )
    
    path("*filtered_10x_mtx/barcodes.tsv.gz", emit: tenx_barcodes)
    path("*filtered_10x_mtx/features.tsv.gz", emit: tenx_features)
    path("*filtered_10x_mtx/matrix.mtx.gz", emit: tenx_matrix)

    path(
      "${outfile}-filtered_10x_mtx-file_list.tsv",
       emit: results_list
    )  
    tuple(
      val(experiment_id),
      val(biopsy_type),
      path("*_unfiltered.h5"),
      path(expected_cells),
      path(total_droplets_include),
      emit: experimentid_outdir_cellbenderunfiltered_expectedcells_totaldropletsinclude
    )
    tuple(
      val(experiment_id),
      val(biopsy_type),
      path(file_10x_barcodes),
      path(file_10x_features),
      path(file_10x_matrix),
      path("*_filtered.h5"),
      emit: experiment_id_cb_plot_input
    )
    tuple(
      val(biopsy_type),
      path(file_10x_barcodes),
      path(file_10x_features),
      path(file_10x_matrix),
      path("*_filtered.h5"),
      emit: cb_plot_input
    )
    
    
    script:
    outfile = "cellbender"
    """
export LD_PRELOAD=/opt/conda/envs/conda_cellbender/lib/libmkl_core.so:/opt/conda/envs/conda_cellbender/lib/libmkl_sequential.so

for i in \$(ls *.h5); do
  echo \$i
  out_file=\$(echo \$i | sed s/"\\."/"pt"/g | sed s/"pth5"/"\\.h5"/)
  # If outfile does not exist, move i to out_file
  [ ! -f \$out_file ] && mv \$i \$out_file
  echo \$out_file
done

python ${projectDir}/../bin/032-clean_cellbender_results.py \\
  --nf_outdir_tag ${params.outdir}/cellbender/3_preprocess_output/${biopsy_type}/${experiment_id} \\
  --cb_outfile_tag ${outfile} \\
  --experiment_id ${experiment_id} \\
  --fpr '${fpr}' \\
  --cb_params ${cb_params}

fprid=\$(echo $fpr | sed s'/\\./pt/'g)
cp -L cellbender_filtered.h5 cellbender_FPR_\${fprid}_filtered.h5
cp -L cellbender_unfiltered.h5 cellbender_FPR_\${fprid}_unfiltered.h5
    """
}
