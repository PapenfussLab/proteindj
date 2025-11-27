process AnalysePredictions {
    label 'pyrosetta_tools'

    publishDir "${params.out_dir}/results", mode: 'copy', pattern: "*.csv"
    publishDir "${params.out_dir}/run/analysis", mode: 'copy', pattern: "analysis.log"

    input:
    path pdbs

    output:
    path "best_designs.jsonl", emit: jsonl, topic: metadata_ch_fold_seq
    path "analysis.log", emit: log

    script:
    def num_processes = task.cpus - 1

    """
    python3 -u /scripts/analyse_best_designs.py \
        --pdb_dir ./ \
        --output "best_designs.jsonl" \
        --verbose \
        --num_processes ${num_processes}
    """
}
process FilterAnalysis {
    label 'python_tools'
    publishDir "${params.out_dir}/run/filter_analysis", mode: 'copy', pattern: '*.log'

    input:
    path(jsonl_file)
    path(pdb_files)

    output:
    path ("output/*.pdb"), emit: pdbs, optional: true
    path "filter_analysis_${task.index}.log"
    path ("filtered.jsonl"), emit: jsonl, optional: true

    script:
    // Format parameters for analysis filtering
    def paramString = Utils.formatFilterParams(
        params,
        "pr",
        [
            "min_helices",
            "max_helices",
            "min_strands",
            "max_strands",
            "min_total_ss",
            "max_total_ss",
            "min_rog",
            "max_rog",
            "min_intface_bsa",
            "min_intface_shpcomp",
            "min_intface_hbonds",
            "max_intface_unsat_hbonds",
            "max_intface_deltag",
            "max_intface_deltagtobsa",
            "min_intface_packstat",
            "max_tem",
            "max_surfhphobics"
        ],
    ) + " " + Utils.formatFilterParams(
        params,
        "seq",
        [
            "min_ext_coef",
            "max_ext_coef",
            "min_pi",
            "max_pi"
        ],
    )

    """    
    python -u /scripts/filter_analysis.py \
        --jsonl-file ${jsonl_file} \
        --pdb-directory ./ \
        ${paramString} \
        --output-directory output \
        2>&1 | tee filter_analysis_${task.index}.log
    """
}
