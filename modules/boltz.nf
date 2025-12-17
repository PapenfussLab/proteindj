process PrepBoltz {
    label 'python_tools'

    input:
    path pdb_files
    path msa_file, stageAs: 'target_msa.a3m'

    output:
    path ("output/*.yaml"), emit: yamls
    path ("output/templates"), emit: templates, optional: true
    path ("output/*.a3m"), emit: msa_file, optional: true

    script:
    // Determine template chain based on design mode
    def templateChain = params.design_mode in ['binder_denovo', 'binder_foldcond', 'binder_motifscaff', 'binder_partialdiff', 'bindcraft_denovo'] ? 'B' : 'A'
    
    // Determine MSA chain based on design mode (same as template chain)
    def msaChain = templateChain
    
    // Build template parameters if enabled
    def templateParams = params.boltz_use_templates ? "--use-template --template-chain ${templateChain}" : ""
    def templateForceParam = params.boltz_use_templates && params.boltz_template_force ? "--template-force" : ""
    def templateThresholdParam = params.boltz_use_templates && params.boltz_template_threshold ? "--template-threshold ${params.boltz_template_threshold}" : ""
    
    // Build MSA parameter if provided
    def msaParam = params.boltz_input_msa ? "--msa-file target_msa.a3m --msa-chain ${msaChain}" : ""
    
    """
    # Generate yaml files containing sequences for Boltz-2 prediction 
    python /scripts/prep_boltz_yaml.py \
        --input "./" \
        --output "output" \
        ${templateParams} \
        ${templateForceParam} \
        ${templateThresholdParam} \
        ${msaParam}
    """
}

process RunBoltz {
    label 'Boltz'
    label 'gpu'
    publishDir "${params.out_dir}/run/boltz", mode: 'copy', pattern: "*.log"
    tag "B${batch_id}"

    input:
    tuple val(batch_id), path(yamls), path(templates_dir, stageAs: 'templates'), path(msa_file)

    output:
    tuple path("predictions/*.pdb"), path("predictions/*.json"), emit: pdbs_jsons   // For AlignBoltz
    path("predictions/*.npz"), emit: npzs
    path ("*.log"), emit: logs

    script:
    """
        # We specify tmp directories as some python packages try to write to the user home directory outside the container
        mkdir tmp
        export NUMBA_CACHE_DIR=tmp
        export XDG_CONFIG_HOME=tmp
        export TRITON_CACHE_DIR=tmp
        export HOME=tmp
        
        mkdir yamls
        for file in \$(find *.yaml); do
            cp -L "\$file" ./yamls/
        done

        boltz predict \
            ./yamls/ \
            --output_format pdb \
            --diffusion_samples ${params.boltz_diffusion_samples} \
            --recycling_steps ${params.boltz_recycling_steps} \
            --sampling_steps ${params.boltz_sampling_steps} \
            ${params.boltz_use_potentials ? '--use_potentials' : ''} \
            --cache /boltzcache \
            ${params.boltz_extra_config ? params.boltz_extra_config : ''} \
            2>&1 | tee boltz_${batch_id}.log
    
        # Move output files out of nested directories and rename to fold_X_seq_X_boltzpred.pdb|json|npz
        mkdir -p predictions
        for dir in boltz_results_yamls/predictions/fold_*_seq_*; do
            # Extract input name from directory path
            inputname=\$(basename "\$dir")
            # Process PDB file
            if [ -f "\${dir}/\${inputname}_model_0.pdb" ]; then
                mv "\${dir}/\${inputname}_model_0.pdb" "predictions/\${inputname}_boltzpred.pdb"
            fi
            # Process JSON file 
            if [ -f "\${dir}/confidence_\${inputname}_model_0.json" ]; then
                mv "\${dir}/confidence_\${inputname}_model_0.json" "predictions/\${inputname}_boltzpred.json"
            fi
            # Process PAE NPZ file
            if [ -f "\${dir}/pae_\${inputname}_model_0.npz" ]; then
                mv "\${dir}/pae_\${inputname}_model_0.npz" "predictions/pae_\${inputname}_boltzpred.npz"
            fi
            # Process PLDDT NPZ file
            if [ -f "\${dir}/plddt_\${inputname}_model_0.npz" ]; then
                mv "\${dir}/plddt_\${inputname}_model_0.npz" "predictions/plddt_\${inputname}_boltzpred.npz"
            fi
        done
    """
}

process AnalyseBoltz {
    label 'python_tools'
    publishDir "${params.out_dir}/run/boltz", mode: 'copy', pattern: "*.log"

    input:
    tuple path(pdb_files), path(json_files)
    path(npz_files)

    output:
    tuple path("output/*.pdb"), path("output/*.json"), emit: pdbs_jsons
    path "analyseboltz_*.log"

    script:
    def num_processes = task.cpus - 1
    """
    # Copy PDB files to output directory (unchanged)
    mkdir -p output
    cp *.pdb output/
    
    # Run interface scoring batch to update JSON files with metrics
    python /scripts/analyse_boltz_batch.py \
        --input-dir ./ \
        --output-dir output \
        --ipsae-script-path /scripts/analyse_boltz_calc.py \
        --pae-cutoff 10 \
        --dist-cutoff 10 \
        --max-workers ${num_processes} \
        --verbose \
        2>&1 | tee analyseboltz_${task.index}.log
    """
}

process AlignBoltz {
    label 'python_tools'
    publishDir "${params.out_dir}/run/boltz", mode: 'copy', pattern: "*.log"

    input:
    tuple path(pdb_files), path(json_files)
    path designs
    val design_type

    output:
    tuple path("aligned/*.pdb"), path("aligned/*.json"), emit: pdbs_jsons
    path "alignment_*.log"
    path ("boltz_metadata_*.jsonl"), topic: metadata_ch_fold_seq

    script:

    def num_processes = task.cpus - 1

    """
    # Script to align predictions to designs and calculate RMSD
    # Also, extracts and renames metadata from json files (includes interface metrics from AnalyseBoltz)
    python /scripts/align_boltz.py \
        --design_dir ./ \
        --boltz_dir ./ \
        --output_dir aligned \
        --design_type ${design_type} \
        --ncpus ${num_processes} \
        2>&1 | tee alignment_${task.index}.log
    
    # metadata convert script to combine json files into jsonl
    python /scripts/metadata_converter.py \
        --converter boltz \
        --input_dir aligned \
        --input_ext .json \
        --output_file boltz_metadata_${task.index}.jsonl
    """
}
process FilterBoltz {
    label 'python_tools'
    publishDir "${params.out_dir}/run/filter_boltz", mode: 'copy', pattern: '*.log'

    input:
    tuple path(pdb_files), path(json_files)

    output:
    path ("output/*.pdb"), emit: pdbs, optional: true
    path ("filter_boltz_${task.index}.log"), emit: log
    path ("filtered.jsonl"), emit: jsonl, optional: true

    script:
    def paramString = Utils.formatFilterParams(
        params,
        "boltz",
        [
            "max_rmsd_overall",
            "max_rmsd_binder",
            "max_rmsd_target",
            "min_conf_score",
            "min_ptm",
            "min_ptm_binder",
            "min_ptm_target",
            "min_ptm_interface",
            "min_plddt",
            "min_plddt_interface",
            "max_pde",
            "max_pde_interface",
            "min_ipSAE_min",
            "min_LIS",
            "min_pDockQ2_min",
            "max_pae_interaction",
        ],
    )

    """
    python -u /scripts/filter_boltz.py \\
        --json-directory ./ \\
        ${paramString} \\
        --output-directory output \\
        2>&1 | tee filter_boltz_${task.index}.log
    """
}
