process PrepBC {
    label 'python_tools'

    input:
    val(batch_id)
    val(batch_size)
    path(input_pdb)
    path(advanced_json)
    
    output:
    tuple path("batch_*.json"), path("advanced_settings.json"), val(batch_id)

    script:
    """
    #!/usr/bin/env python3
    
    import os
    import json
    from Bio.PDB import PDBParser
    
    starting_pdb = os.path.basename("${input_pdb}")
    design_path = "designs"
    binder_name = "batch_${batch_id}"

    # Function to extract protein chains from PDB
    def get_protein_chains(pdb_file):
        # Extract protein chain IDs from PDB file.
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        
        protein_chains = []
        for model in structure:
            for chain in model:
                # Check if chain has at least one standard amino acid residue
                has_protein = False
                for residue in chain:
                    if residue.id[0] == ' ':  # Standard amino acid (not hetero)
                        has_protein = True
                        break
                
                if has_protein and chain.id not in protein_chains:
                    protein_chains.append(chain.id)
        
        return sorted(protein_chains)

    # Get chains - either from parameter or auto-detect
    # Check for both empty string and literal "null" from Nextflow
    bc_chains_param = "${params.bc_chains}".strip()
    if bc_chains_param and bc_chains_param != "null":
        # Use provided chains
        chains_list = [c.strip() for c in "${params.bc_chains}".split(',')]
        print(f"Using provided chains: {chains_list}")
    else:
        # Auto-detect chains from PDB
        chains_list = get_protein_chains("${input_pdb}")
        print(f"Auto-detected chains from PDB: {chains_list}")
        if not chains_list:
            raise ValueError("No protein chains found in PDB file. Please provide chains explicitly with bc_chains parameter.")
    
    # Convert to comma-separated string for BindCraft
    chains = ",".join(chains_list)
    
    # Parse hotspots - convert to comma-separated string for BindCraft
    hotspots_param = "${params.hotspot_residues}".strip()
    if hotspots_param and hotspots_param != "null":
        hotspots = ",".join([h.strip() for h in hotspots_param.split(',') if h.strip()])
    else:
        hotspots = ""
    
    # Parse lengths - handle single value or range, store as list for JSON array
    length_str = "${params.design_length}".strip()
    if '-' in length_str:
        lengths = [int(x.strip()) for x in length_str.split('-')]
    else:
        lengths = [int(length_str)]

    # Create settings dictionary
    # BindCraft expects chains and hotspots as comma-separated strings
    # but lengths as an array
    settings = {
        "design_path": design_path,
        "binder_name": binder_name,
        "starting_pdb": starting_pdb,
        "chains": chains,  # String: "A,B"
        "target_hotspot_residues": hotspots,  # String: "A50,B75"
        "lengths": lengths,  # List: [60, 80]
        "number_of_final_designs": ${batch_size}
    }

    # Write settings JSON file
    target_settings_path = os.path.join(binder_name + ".json")
    with open(target_settings_path, 'w') as f:
        json.dump(settings, f, indent=4)
    print(f"Wrote settings JSON file: {target_settings_path}")
    print(f"  Chains: {chains}")
    print(f"  Hotspots: {hotspots}")
    print(f"  Lengths: {lengths}")

    # Load and modify advanced settings JSON
    advanced_json_path = "${advanced_json}"
    with open(advanced_json_path, 'r') as f:
        advanced_settings = json.load(f)
    
    # Update amino acid types to avoid during design
    omit_AAs_param = "${params.bc_omit_AAs}".strip()
    if omit_AAs_param and omit_AAs_param != "null":
        omit_AAs = ",".join([aa.strip() for aa in omit_AAs_param.split(',') if aa.strip()])
    else:
        omit_AAs = ""
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
    print(f"Modified advanced settings JSON file: {modified_advanced_path}")
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
    label 'python_tools'
    publishDir "${params.out_dir}/run/bc", mode: 'copy', pattern: "analysis_*.log"

    input:
    path(pdbs_csvs)

    output:
    tuple path("processed/*.pdb"), path("processed/*.json"), emit: pdbs_jsons
    path "bindcraft_analysis.log"
    path("bindcraft_metadata.jsonl"), topic: metadata_ch_fold

    script:
    """
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