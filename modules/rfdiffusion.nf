process GenerateRFDContigs {
    label 'python_tools'

    input:
    path(input_pdb)
    val(design_mode)
    val(is_binder_mode)
    
    output:
    env 'CONTIGS'
    
    script:
    """
    CONTIGS=\$(python /scripts/generate_contigs.py ${input_pdb} \
        ${params.design_length ? "--design_length ${params.design_length}" : ''} \
        --design_mode ${design_mode} \
        ${is_binder_mode ? '--is_binder' : ''})
    echo "Generated contigs: \$CONTIGS"
    export CONTIGS
    """
}
process GenerateRFDFoldCond {
    label 'python_tools'

    input:
    path(input_pdb)
    
    output:
    path("target_ss.pt"), emit: target_ss
    path("target_adj.pt"), emit: target_adj
    
    script:
    """
    python /scripts/create_scaffolds.py --input_pdb ${input_pdb} --out_dir "./" --num_processes 1
    mv *_ss.pt target_ss.pt
    mv *_adj.pt target_adj.pt
    """
}
process RunRFD {
    label 'RFDiffusion'
    label 'gpu'
    tag "B${batch_id}"

    publishDir "${params.out_dir}/run/rfd", mode: 'copy', pattern: "*.log"

    beforeScript """
        mkdir -p outputs schedules .dgl
    """

    input:
    tuple val(batch_id), val(batch_size), val(design_startnum), val(mode), path(input_files), val(rfdCommand), path(target_adj), path(target_ss)

    output:
    path ("rfd_results/*.pdb"), emit: pdbs
    tuple path("rfd_results/*.pdb"), path("rfd_results/*.json"), emit: pdbs_jsons
    path "*.log"
    path ("rfd_metadata_${batch_id}.jsonl"), topic: metadata_ch_fold

    script:
    // Note: rfdCommand is precomputed in the workflow body (main.nf) rather than here, so
    // that the task cache is keyed only on this specific string rather than the entire
    // global params object.
    def inference_log_filename = "rfd_${batch_id}.log"
    
    """
    echo "Running RFdiffusion for batch ${batch_id} in ${mode} mode"
    echo "RFdiffusion command: ${rfdCommand} inference.num_designs=${batch_size}"
    python3 ${rfdCommand} \
        inference.model_directory_path=/app/RFdiffusion/models \
        inference.schedule_directory_path=/app/RFdiffusion/schedules \
        inference.design_startnum=${design_startnum} \
        inference.num_designs=${batch_size} 2>&1 | tee ${inference_log_filename}
    
    python3 /scripts/metadata_converter.py --input_dir rfd_results --converter rfd --input_ext trb -o rfd_metadata_${batch_id}.jsonl
    """
}
process FilterFold {
    label 'python_tools'
    publishDir "${params.out_dir}/run/filter_fold", mode: 'copy', pattern: "*.log"

    input:
    tuple path(pdb_files), path(json_files)
    val(paramString)

    output:
    tuple path("filtered_output/*.pdb"), path("filtered_output/*.json"), emit: pdbs_jsons, optional: true
    path ("fold_data_*.jsonl"), topic: metadata_ch_fold
    path "filter_fold_*.log"

    script:

    def num_processes = task.cpus - 1
    """
    python /scripts/filter_fold.py \
        --input-dir . \
        --output-dir "filtered_output" \
        ${paramString} \
        --ncpus ${num_processes}
    """
}
