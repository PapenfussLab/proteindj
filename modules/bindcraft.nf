process PrepBC {
    cpus 1
    memory '1 GB'

    input:
    val(batch_id)
    val(batch_size)
    path(advanced_json)
    
    output:
    tuple path("batch_*.json"), path("advanced_settings.json"), val(batch_id)

    script:
    """
    #!/usr/bin/env python3
    
    import os
    import json
    
    starting_pdb = os.path.basename("${params.input_pdb}")
    design_path = "designs"
    binder_name = "batch_${batch_id}"

    # Properly convert comma-separated parameters to Python lists
    chains = ",".join([c.strip() for c in "${params.bc_chains}".split(',')]) if "${params.bc_chains}" else ""
    hotspots = ",".join([h.strip() for h in "${params.hotspot_residues}".split(',')]) if "${params.hotspot_residues}" else ""
    lengths = [int(x.strip()) for x in "${params.bc_design_length}".split(',')]

    # Number of final designs is irrelevant as we skip MPNN-AF2 validation
    settings = {
        "design_path": design_path,
        "binder_name": binder_name,
        "starting_pdb": starting_pdb,
        "chains": chains,
        "target_hotspot_residues": hotspots,
        "lengths": lengths,
        "number_of_final_designs": ${batch_size}
    }

    target_settings_path = os.path.join(binder_name+".json")
    with open(target_settings_path, 'w') as f:
        json.dump(settings, f, indent=4)
    print("Wrote settings JSON file")

    # Load and modify advanced settings JSON
    advanced_json_path = "${advanced_json}"
    with open(advanced_json_path, 'r') as f:
        advanced_settings = json.load(f)
    
    # Update amino acid types to avoid during design
    omit_AAs = ",".join([aa.strip() for aa in "${params.bc_omit_AAs}".split(',')]) if "${params.bc_omit_AAs}" else ""
    advanced_settings["omit_AAs"] = omit_AAs
    # Disable BindCraft's internal MPNN routine as we will perform this in ProteinDJ
    advanced_settings["enable_mpnn"] = False
    # Set max trajectories to batch size. BindCraft will run until it achieves this many designs.
    advanced_settings["max_trajectories"] = ${batch_size}
    # Set path to AF2 params directory to bind location in container
    advanced_settings["af_params_dir"] = "/af2params"
    
    # Write the modified advanced settings
    modified_advanced_path = "advanced_settings.json"
    with open(modified_advanced_path, 'w') as f:
        json.dump(advanced_settings, f, indent=4)
    print("Modified advanced settings JSON file")
    """
}

process RunBC {
    label 'BC'
    label 'gpu'

    publishDir "${params.out_dir}/run/bc", mode: 'copy', pattern: "*.log"
    publishDir "${params.out_dir}/run/bc", mode: 'copy', pattern: "designs/trajectory_stats.csv"

    input:
    tuple path(settings_json), path(advanced_json), val(batch_id)
    path(filters_json)
    path(input_pdb)

    output:
    path "*.log"
    tuple path("designs/Trajectory/Relaxed/*.pdb"), path("batch_*.csv"), emit: pdbs_csvs, optional: true

    script:
    """
    # We specify a tmp directory as some python packages try to write to the user home directory outside the container
    mkdir tmp
    export XDG_CACHE_HOME="tmp"
    export MPLCONFIGDIR="tmp"
    
    python /opt/BindCraft/bindcraft.py \
        --settings ${settings_json} \
        --filters ${filters_json} \
        --advanced ${advanced_json} \
        2>&1 | tee bindcraft_${task.index}.log
    
    # Rename CSV file if pdb file exists
    # TODO: Consider the case of 100% rejection and if batch fails completely
    pdb_file=\$(ls designs/Trajectory/Relaxed/*.pdb 2>/dev/null)
    if [[ -n "\${pdb_file}" ]]; then
        mv designs/trajectory_stats.csv batch_${batch_id}.csv
    else
        echo "All designs were rejected"
    fi

    """
}
process AnalyseBC {
    label 'pyrosetta_tools'
    publishDir "${params.out_dir}/run/bc", mode: 'copy', pattern: "analysis_*.log"

    input:
    path(pdbs_csvs)

    output:
    tuple path("processed/*.pdb"), path("processed/*.json"), emit: pdbs_jsons
    path "bindcraft_analysis.log"
    path("bindcraft_metadata.jsonl"), topic: metadata_ch_fold

    script:
    """
    export MAMBA_ROOT_PREFIX=/opt/conda/
    
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate pyrosetta

    # Script to process BindCraft output PDBs and swap chains: A→B, B→A. Will create a new unique fold_id.
    # Also, extracts and renames metadata from CSV to json files
    python /scripts/analyse_bindcraft.py \
        --input_dir ./ \
        --output_dir processed \
        ${params.bc_fix_interface_residues ? '--fix_interface_residues' : ''} \
        2>&1 | tee bindcraft_analysis.log
    
    python /scripts/metadata_converter.py \
        --converter bc \
        --input_dir processed \
        --input_ext .json \
        --output_file bindcraft_metadata.jsonl
    """
}