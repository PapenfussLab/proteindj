process PrepMPNN {
    label 'python_tools'
    publishDir "${params.out_dir}/run/mpnn", mode: 'copy', pattern: "mpnn_prep_*.log"

    input:
    tuple path(pdb_files), path(json_files)

    output:
    path ("mpnn_input/*.pdb"), emit: pdbs
    path ("mpnn_prep_*.log")

    script:
    """
    python /scripts/prep_mpnn_designs.py \
        --input_dir "./" \
        --out_dir "mpnn_input"
    
    # Add unique ID to mpnn prep logfile
    cp mpnn_prep.log mpnn_prep_${task.index}.log

    """
}

process RunMPNN {
    label 'MPNN'
    cpus 2

    publishDir "${params.out_dir}/run/mpnn", mode: 'copy', pattern: "*.log"

    // retry up to 3 times because proteinmpn occasionally fails due to memory leaks
    maxRetries 3

    input:
    tuple val(batch_id), path(pdbs)

    output:
    tuple path("results/*.pdb"), path("results/*.json"), emit: pdbs_jsons
    path ("mpnn_metadata_${batch_id}.jsonl"), topic: metadata_ch_fold_seq
    path "*.log"

    script:
    """
    export OPENBLAS_NUM_THREADS=1 
    export MKL_NUM_THREADS=1

    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate mpnn
    mkdir results
        
    python /dl_binder_design/mpnn_fr/dl_interface_design_multi.py \
        -pdbdir "./" \
        -outpdbdir "./results" \
        -augment_eps ${params.mpnn_backbone_noise} \
        -checkpoint_path "/dl_binder_design/mpnn_fr/ProteinMPNN/${params.mpnn_checkpoint_type}_model_weights/${params.mpnn_checkpoint_model}.pt" \
        -omit_AAs ${params.mpnn_omitAAs} \
        -relax_max_cycles ${params.mpnn_relax_max_cycles} \
        ${params.mpnn_relax_output ? '-relax_output' : ''} \
        -relax_seqs_per_cycle  ${params.mpnn_relax_seqs_per_cycle} \
        -relax_convergence_rmsd ${params.mpnn_relax_convergence_rmsd} \
        -relax_convergence_score ${params.mpnn_relax_convergence_score} \
        -relax_convergence_max_cycles ${params.mpnn_relax_convergence_max_cycles}\
        -seqs_per_struct ${params.seqs_per_design} \
        -temperature ${params.mpnn_temperature} \
        -debug \
        ${params.mpnn_extra_config ? params.mpnn_extra_config : ''} \
        2>&1 | tee mpnn_${batch_id}.log

    python /scripts/metadata_converter.py --input_dir results --input_ext ".json" \
        --converter mpnn  --output_file "mpnn_metadata_${batch_id}.jsonl"
    """
}
process FilterMPNN {
    label 'python_tools'

    publishDir "${params.out_dir}/run/filter_mpnn", mode: 'copy', pattern: '*.log'

    input:
    tuple path(pdb_files), path(json_files)

    output:
    path ("filtered_output/*.pdb"), emit: pdbs, optional: true
    path ("filtered_output/*.json"), emit: jsons, optional: true
    path ("filter_mpnn_${task.index}.log"), emit: logs

    script:
    // Only pass parameters if filter values are provided
    def mpnnParam = Utils.formatFilterParams(params, "mpnn", ["max_score"])

    """    
    python /scripts/filter_mpnn.py \
        --jsons ./ \
        --pdbs ./ \
        ${mpnnParam} \
        --output-dir filtered_output \
        2>&1 | tee filter_mpnn_${task.index}.log
    """
}
