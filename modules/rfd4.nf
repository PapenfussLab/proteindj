process PrepRFD4 {
    label 'python_tools'

    input:
    path(input_pdb)

    output:
    tuple path('prepared/binder.json'), path('prepared/target.cif')

    script:
    """
    python /scripts/prepare_rfd4.py --input-pdb ${input_pdb} \
        --length '${params.design_length}' \
        --hotspots '${params.hotspot_residues ?: ''}' \
        --flexible-residues '${params.flexible_residues ?: ''}' \
        --output-dir prepared
    """
}

process RunRFD4 {
    label 'RFD4'
    label 'gpu'
    tag "B${batch_id}"
    publishDir "${params.out_dir}/run/rfd4", mode: 'copy', pattern: '*.log'

    input:
    tuple val(batch_id), val(designs_per_gpu), path(spec), path(target)

    output:
    tuple val(batch_id), val(designs_per_gpu), path('rfd4_run/inference_samples'), path(target), emit: samples
    path '*.log'

    script:
    // num_replicates is the pipeline's designs per GPU. RFD4's own batch_size
    // controls concurrent inference within this job and cannot exceed that count.
    def inference_batch_size = Math.min(Integer.parseInt(params.rfd4_inference_batch_size.toString()), designs_per_gpu.intValue())
    def triton_cache_dir = params.rfd4_triton_cache ? '/rfd4_triton_cache' : '.rfd4_triton_cache'
    """
    set -o pipefail
    mkdir -p rfd4_run ${triton_cache_dir}/cuequivariance ${triton_cache_dir}/nvrtc ${triton_cache_dir}/triton
    export CUEQ_TRITON_CACHE_DIR=${triton_cache_dir}/cuequivariance
    export CUEQUIVARIANCE_OPS_NVRTC_CACHE_DIR=${triton_cache_dir}/nvrtc
    export TRITON_CACHE_DIR=${triton_cache_dir}/triton
    rfd4-generate hf_ckpt_path=/rfd4_models inputs=${spec} num_replicates=${designs_per_gpu} \
        batch_size=${inference_batch_size} out_dir=rfd4_run \
        compile_lengths=null '~metrics' 2>&1 | tee rfd4_${batch_id}.log
    """
}

process AnalyseRFD4 {
    label 'python_tools'

    input:
    tuple val(batch_id), val(designs_per_gpu), path(samples), path(target)

    output:
    tuple path('processed/fold_*.pdb'), path('processed/fold_*.json'), emit: pdbs_jsons
    path('processed/rfd4_metadata.jsonl'), topic: metadata_ch_fold

    script:
    def first_id = batch_id * designs_per_gpu
    """
    python /scripts/analyse_rfd4.py --samples-dir ${samples} --target-cif ${target} \
        --start-id ${first_id} --output-dir processed
    """
}

process UseRFD4Sequences {
    label 'python_tools'

    input:
    tuple path(pdbs), path(jsons)

    output:
    path('sequences/fold_*_seq_0.pdb'), emit: pdbs, optional: true
    path('sequences/rfd4_seq_metadata.jsonl'), topic: metadata_ch_fold_seq

    script:
    def limits = [
        ['min-ext', params.seq_min_ext_coef], ['max-ext', params.seq_max_ext_coef],
        ['min-pi', params.seq_min_pi], ['max-pi', params.seq_max_pi]
    ].findAll { it[1] != null }.collect { "--${it[0]} ${it[1]}" }.join(' ')
    """
    python /scripts/use_rfd4_sequences.py --input-dir . --output-dir sequences ${limits}
    """
}
