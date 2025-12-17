process RunAF2 {
    label 'AF2'
    label 'gpu'

    publishDir "${params.out_dir}/run/af2", mode: 'copy', pattern: "*.log"

    tag "B${batch_id}"

    input:
    tuple val(batch_id), path(pdbs)

    output:
    tuple path("outputs/*.pdb"), path("*.json"), emit: pdbs_jsons
    path ("*.json"), emit: json, topic: metadata_ch_fold_seq
    path "*.log"

    script:
    """
        python3 /dl_binder_design/af2_initial_guess/predict.py \
            -pdbdir ./ \
            -outpdbdir outputs \
            -scorefilename ${batch_id}_out.sc \
            ${params.af2_initial_guess ? '' : '-no_initial_guess '} \
            ${params.af2_extra_config ? params.af2_extra_config : ''} \
            2>&1 | tee af2_${batch_id}.log
        
        python3 /scripts/metadata_converter.py \
            --input_files ${batch_id}_out.sc \
            --converter af2 \
            --split_by_description
    """
}
process AlignAF2 {
    label 'python_tools'
    publishDir "${params.out_dir}/run/align", mode: 'copy', pattern: "alignment_*.log"

    input:
    path af2_pdbs
    path reference_pdb

    output:
    path "aligned/*.pdb", emit: pdbs, optional: true
    path "alignment_${task.index}.log", emit: logs, optional: true

    script:

    def num_processes = task.cpus - 1

    """
    mkdir -p referencepdb
    mv ${reference_pdb} referencepdb/.

    python /scripts/align_af2.py \
        --input_dir ./ \
        --output_dir "aligned" \
        --reference referencepdb/*.pdb \
        --ncpus ${num_processes} \
        2>&1 | tee alignment_${task.index}.log
            
    """
}
process FilterAF2 {
    label 'python_tools'
    publishDir "${params.out_dir}/run/filter_af2", mode: 'copy', pattern: '*.log'

    input:
    tuple path(pdb_files), path(json_files)

    output:
    path ("output/*.pdb"), emit: pdbs, optional: true
    path "filter_af2_${task.index}.log"
    path ("filtered.jsonl"), emit: jsonl, optional: true

    script:
    // Only pass parameters if filter values are provided
    def paramString = Utils.formatFilterParams(
        params,
        "af2",
        [
            "max_pae_interaction",
            "max_pae_overall",
            "max_pae_binder",
            "max_pae_target",
            "min_plddt_overall",
            "min_plddt_binder",
            "min_plddt_target",
            "max_rmsd_overall",
            "max_rmsd_binder_bndaln",
            "max_rmsd_binder_tgtaln",
            "max_rmsd_target"
        ],
    )

    """    
    python -u /scripts/filter_af2.py \
        --json-directory ./ \
        ${paramString} \
        --output-directory output \
        2>&1 | tee filter_af2_${task.index}.log
    """
}
