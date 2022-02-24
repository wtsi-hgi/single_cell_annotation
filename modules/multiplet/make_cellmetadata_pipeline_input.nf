process make_cellmetadata_pipeline_input {
    // Makes a input tsv file for the main pipeline.
    // ------------------------------------------------------------------------
    
    publishDir  path: "${params.outdir}/multiplet/2_make_cellmetadata_pipeline_input",
        mode: "${params.copy_mode}",
        overwrite: "true"
    
    input:
    path("*multiplet_calls_published.txt")
    
    output:
    path('file_cellmetadata.tsv', emit: file__cellmetadata)
    
    script:
"""
# Note: the default paste delim is tab
cat *multiplet_calls_published.txt \
    | awk 'BEGIN{print "experiment_id\tdata_path_cellmetadata"}1' \
    > file_cellmetadata.tsv
"""
}
