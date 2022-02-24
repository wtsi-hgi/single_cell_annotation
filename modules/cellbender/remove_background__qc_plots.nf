process remove_background__qc_plots {

  publishDir  path: "${params.outdir}/cellbender/4.1_qc_plots/${biopsy_type}/${experiment_id}",
      saveAs: {filename ->
        if (filename.equalsIgnoreCase("barcodes.tsv.gz")) {
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
	val(biopsy_type),
	path(file_10x_barcodes),
	path(file_10x_features),
	path(file_10x_matrix),
	file(h5_filtered_cellbender)
    )

  output:
    path("plots/*.png") optional true
    path("plots/*.pdf") optional true

  script:
    h5_filtered_cellbender = h5_filtered_cellbender.join(",")
    """
    mkdir -p txd_input
    ln --physical ${file_10x_barcodes} txd_input/barcodes.tsv.gz
    ln --physical ${file_10x_features} txd_input/features.tsv.gz
    ln --physical ${file_10x_matrix} txd_input/matrix.mtx.gz

    # Make a file with list of our files
    echo "${h5_filtered_cellbender}" | sed s/","/"\\n"/g > files.txt

    for i in \$(cat files.txt); do
    echo \$i
    out_file=\$(echo \$i | sed s/".h5"//)
    python ${projectDir}/../bin/035-analyse_cellbender_results.py \\
      --tenxdata_path txd_input  \\
      --h5_cellbender \$i \\
      --output_file cellbender_results-\$out_file \\
      --number_cpu ${task.cpus}
    done

    rm files.txt
    mkdir -p plots
    mv *pdf plots/ 2>/dev/null || true
    mv *png plots/ 2>/dev/null || true
    """
}
