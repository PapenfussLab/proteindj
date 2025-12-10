process PrepFAMPNN {
    label 'pyrosetta_tools'

    input:
    tuple path(pdb_files), path(json_files)

    output:
    path ("fampnn_input/*.pdb"), emit: pdbs
    path ("*.csv"), emit: csv

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate pyrosetta

    # Restore missing side-chains required by FAMPNN    
    python /scripts/prep_fampnn_designs.py \
        --input_dir "./" \
        --out_dir "fampnn_input"
    
    # Define residue ranges to be fixed during sequence generation in CSV
    python /scripts/prep_fampnn_csv.py \
        --input_dir "./" \
        --out_csv "fampnn.csv" \
        ${params.fampnn_fix_target_sidechains ? "--fix_target_sidechains" : ''}
    """
}

process RunFAMPNN {
    label 'FAMPNN'
    label 'gpu'

    publishDir "${params.out_dir}/run/fampnn", mode: 'copy', pattern: "*.log"

    input:
    tuple val(batch_id), path(pdbs), path(csv)
    val analysis_chain_id

    output:
    tuple path("results/*.pdb"), path("results/*.json"), emit: pdbs_jsons
    path ("fampnn_metadata_${batch_id}.jsonl"), topic: metadata_ch_fold_seq
    path "*.log"

    script:
    """
    mkdir -p results

    python /app/fampnn/fampnn/inference/seq_design.py \
        batch_size=16 \
        checkpoint_path=/app/fampnn/weights/fampnn_0_3.pt \
        exclude_cys=${params.fampnn_exclude_cys} \
        fixed_pos_csv=${csv} \
        num_seqs_per_pdb=${params.seqs_per_design} \
        pdb_dir="./" \
        presort_by_length=true \
        psce_threshold=${params.fampnn_psce_threshold}  \
        temperature=${params.fampnn_temperature} \
        out_dir="fampnn_output" \
        ${params.fampnn_extra_config ? params.fampnn_extra_config : ''} \
        2>&1 | tee fampnn_${task.index}.log

    # Rename output files from fold_X_sampleY.pdb to fold_X_seq_Y.pdb
    for file in fampnn_output/samples/fold_*_sample*.pdb; do
        # Extract the base filename
        base_name=\$(basename "\$file")
        new_name=\$(echo "\$base_name" | sed 's/sample/seq_/')
        cp "\$file" "results/\$new_name"
    done

    python /scripts/analyse_fampnn.py \\
        --input_dir results \
        --chain_id ${analysis_chain_id} \
        --ignore_cbeta \
        --out_dir results

    # Combine metadata to jsonl file
    python /scripts/metadata_converter.py --input_dir results --input_ext ".json" \
        --converter fampnn --output_file "fampnn_metadata_${batch_id}.jsonl"
    
    """
}
process FilterFAMPNN {
    label 'python_tools'

    publishDir "${params.out_dir}/run/filter_fampnn", mode: 'copy', pattern: '*.log'

    input:
    tuple path(pdb_files), path(json_files)

    output:
    path ("filtered_output/*.pdb"), emit: pdbs, optional: true
    path ("filtered_output/*.json"), emit: jsons, optional: true
    path ("filter_fampnn_${task.index}.log"), emit: logs

    script:
    // Only pass parameters if filter values are provided
    def fampnnParam = Utils.formatFilterParams(params, "fampnn", ["max_psce"])

    """    
    python /scripts/filter_fampnn.py \
        --jsons ./ \
        --pdbs ./ \
        ${fampnnParam} \
        --output-dir filtered_output \
        2>&1 | tee filter_fampnn_${task.index}.log
    """
}