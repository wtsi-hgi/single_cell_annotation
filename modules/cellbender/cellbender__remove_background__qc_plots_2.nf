

process cellbender__remove_background__qc_plots_2 {


  // Second set of QC plots from cellbdender.
  // This task compare the cellbender output with both the cellranger filtered and cellragner raw outputs
  // ------------------------------------------------------------------------
  tag { "$experiment_id" }
  publishDir "${params.outdir}/cellbender/4.2_qc_plots/vs_cellranger/${biopsy_type}/$experiment_id", pattern: "fpr_${fpr}/${experiment_id}/*.png", 
    saveAs: {filename ->
    filename.replaceAll("fpr_${fpr}/${experiment_id}/", "fpr_${fpr}/")
  },
    mode: "${params.copy_mode}",
    overwrite: "true"
  
  input:
    tuple val(experiment_id),
	path(cellbender_unfiltered_h5s),
	path(expectedcells),
	path(totaldropletsinclude),
	path(raw_cellranger_mtx),
	path(filtered_cellranger_mtx),
	val(fpr),
	val(biopsy_type)
    // tuple val(experiment_id), val(outdir), path(cellbender_unfiltered_h5s), path(expectedcells), path(totaldropletsinclude), path(raw_cellranger_mtx), path(filtered_cellranger_mtx), val(fpr)

    output:
    path("fpr_${fpr}/${experiment_id}/*.png"), emit: plots_png 

  script:
  """
  fprid=\$(echo $fpr | sed s'/\\./pt/'g)
  cellbender_unfiltered_h5=cellbender_FPR_\${fprid}_unfiltered.h5

  n_expected_cells=\$(cat $expectedcells)
  n_total_droplets_included=\$(cat $totaldropletsinclude)

  python ${projectDir}/../bin/037-plot_cellranger_vs_cellbender.py \\
    --samplename \"${experiment_id}\" \\
    --raw_cellranger_mtx \"${raw_cellranger_mtx}\" \\
    --filtered_cellranger_mtx \"${filtered_cellranger_mtx}\" \\
    --cellbender_unfiltered_h5 \"\$cellbender_unfiltered_h5\" \\
    --fpr \"${fpr}\" \\
    --n_expected_cells \"\${n_expected_cells}\" \\
    --n_total_droplets_included \"\${n_total_droplets_included}\" \\
    --out_dir \$PWD
  """
}
