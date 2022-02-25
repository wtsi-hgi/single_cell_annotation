process merge {

    tag { "${biopsy_type}" }
    publishDir  path: "${params.outdir}/merge/1_merge/${biopsy_type}/",
        mode: "${params.copy_mode}",
        overwrite: "true"
    
    input:
    tuple(
        val(biopsy_type),
	path(input_anndata_files))
    path(samples_metainfo_tsv)
    path(filter_params)
    
    output:
    tuple(
        val(biopsy_type),
        path("*.h5ad"),
        emit: biopsy_type_h5ad)
    tuple(
        val(biopsy_type),
        path("*-cell_filtered_per_experiment.tsv.gz"),
        emit: cell_filtered_per_experiment_tsv_gz)
    path("input_file_paths.tsv")
    
    script:
"""
ls . | grep 'h5ad\$' | sort > h5ad_paths.txt
cat h5ad_paths.txt | sed s'/___.*\$//'g > sample_names.txt
printf \"experiment_id\\tanndata_file\\n" > input_file_paths.tsv
paste -d \"\\t\" sample_names.txt h5ad_paths.txt >> input_file_paths.tsv

Rscript ${projectDir}/../bin/041-clean_metadata.R \\
  --input_file ${samples_metainfo_tsv}

python ${projectDir}/../bin/042-scanpy_merge.py \\
    --number_cpu ${task.cpus} \\
    --anndata_files input_file_paths.tsv \\
    --params_yaml ${filter_params} \\
    --sample_metadata_file samples_metainfo-cleaned.tsv.gz \\
    --metadata_key ${params.merge.metadata_key} \\
    --anndata_compression_opts ${params.merge.anndata_compression_opts}
"""
}
