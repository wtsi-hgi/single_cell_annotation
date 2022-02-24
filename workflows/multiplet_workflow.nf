include { run_scrublet } from '../modules/multiplet/run_scrublet.nf'
include { make_cellmetadata_pipeline_input } from '../modules/multiplet/make_cellmetadata_pipeline_input.nf'

workflow multiplet_workflow {
    take:
        channel__file_paths_10x_scrublet_params
        //expected_multiplet_rate
        //n_simulated_multiplet
        //multiplet_threshold_method
        //scale_log10
    main:
    // Identify multiplets using scrublet.
    run_scrublet(
        channel__file_paths_10x_scrublet_params)
    
    make_cellmetadata_pipeline_input(
        run_scrublet.out.multiplet_calls_published.collect())
    
    emit:
    keras_input = run_scrublet.out.keras_input
    // Return merged input data file.
    file__cellmetadata = make_cellmetadata_pipeline_input.out.file__cellmetadata
    multiplet_calls = run_scrublet.out.multiplet_calls
    
}
