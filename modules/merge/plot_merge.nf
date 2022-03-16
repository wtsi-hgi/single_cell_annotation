process plot_merge {

    tag { "${biopsy_type}" }
    publishDir  path: "${params.outdir}/merge/2_plot/${biopsy_type}/",
        mode: "${params.copy_mode}",
        overwrite: "true"
    
    input:
    tuple(
        val(biopsy_type),
        path(merged_h5ad),
        path(tsv_file_gz))
    
    output:
    tuple(
        val(biopsy_type),
        path("*.png"),
        path("*predicted_sex_not_match_annotated.tsv"),
        emit: out)
    
    script:
"""
python ${projectDir}/../bin/043-plot_filtered_cells.py \\
  --tsv_file ${tsv_file_gz}

python ${projectDir}/../bin/043-plot_final_data.py \\
  --h5ad_file ${merged_h5ad}

# plot to check for sample swaps:
python ${projectDir}/../bin/028-plot_sex.py \\
  --h5_anndata ${merged_h5ad}
"""
}
