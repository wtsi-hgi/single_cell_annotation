process plot_merge {

    tag { "${biopsy_type}" }
    publishDir  path: "${params.outdir}/merge/2_plot/${biopsy_type}/",
        mode: "${params.copy_mode}",
        overwrite: "true"
    
    input:
    tuple(
        val(biopsy_type),
        path(tsv_file_gz))
    
    output:
    tuple(
        val(biopsy_type),
        path("*.png"),
        emit: out)
    
    script:
"""
python ${projectDir}/../bin/043-plot_filtered_cells.py \\
  --tsv_file ${tsv_file_gz}
"""
}

//--tsv_file /lustre/scratch119/realdata/mdt2/projects/sc-eqtl-ibd/analysis/leland_dev/cdi/adata-cell_filtered_per_experiment.tsv.gz
