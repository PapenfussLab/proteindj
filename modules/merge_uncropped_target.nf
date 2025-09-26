process MergeUncroppedTarget {
    label 'pyrosetta_tools'

    input:
    path pdb_files
    path uncropped_target_pdb

    output:
    path ("merged_pdbs/*.pdb"), emit: pdbs

    script:

    def num_processes = task.cpus - 1

    """
    mkdir -p uncropped_pdb
    mv ${uncropped_target_pdb} uncropped_pdb/. 

    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate pyrosetta

    python /scripts/merge_uncropped_target.py \
    --uncropped_pdb uncropped_pdb/*.pdb \
    --input_dir './'  \
    --output_dir merged_pdbs \
    --ncpus ${num_processes}
    """
}