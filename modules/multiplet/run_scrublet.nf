process run_scrublet {
    // Runs scrublet for each sample.
    // ------------------------------------------------------------------------
    //tag { output_dir }

    publishDir  path: "${params.outdir}/multiplet/1_scrublet/${biopsy_type}/${experiment_id}",
                saveAs: {filename ->
                    if (filename.endsWith("multiplet_calls_published.txt")) {
                        null
                    } else {
                        filename.replaceAll("", "")
                    }
                },
                mode: "${params.copy_mode}",
                overwrite: "true"

    //each smplid_and_datapath
    input:
    tuple(
        val(experiment_id),
	val(biopsy_type),
        path(file_10x_barcodes),
        path(file_10x_features),
        path(file_10x_matrix),
        val(expected_multiplet_rate),
        val(n_simulated_multiplet),
        val(multiplet_threshold_method),
        val(scale_log10))

    output:
    tuple(
        val(experiment_id),
	val(biopsy_type),
        path("*-scrublet.h5ad"),
        emit: keras_input)

    path("${experiment_id}-scrublet.tsv.gz", emit: multiplet_calls)
    path(
        "${experiment_id}-multiplet_calls_published.txt",
        emit: multiplet_calls_published
    )
    path("plots/*.pdf") optional true
    path("plots/*.png") optional true

    script:
        // Check to see if we should use use log10 of the doublet simulations
        // to derive the threshold
        cmd__scale_log10 = ""
        if (scale_log10 == "True") {
            cmd__scale_log10 = "--scale_log10"
        }
        """
        rm -fr plots
        TMP_DIR=\$(mktemp -d -p \$(pwd))
        ln --physical ${file_10x_barcodes} \$TMP_DIR
        ln --physical ${file_10x_features} \$TMP_DIR
        ln --physical ${file_10x_matrix} \$TMP_DIR
        python ${projectDir}/../bin/0015-run_scrublet-cdi_version.py \
            --tenxdata_dir \$TMP_DIR \
            --expected_multiplet_rate ${expected_multiplet_rate} \
            --n_simulated_multiplet ${n_simulated_multiplet} \
            --multiplet_threshold_method ${multiplet_threshold_method} \
            ${cmd__scale_log10} \
            --output_file ${experiment_id}
        echo -e "${experiment_id}\t${params.outdir}/multiplet/1_scrublet/${experiment_id}/${experiment_id}-scrublet.tsv.gz" > \
            ${experiment_id}-multiplet_calls_published.txt
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}
