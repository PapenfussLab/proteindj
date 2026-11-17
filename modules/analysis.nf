process AnalysePredictions {
    label 'python_tools'

    publishDir "${params.out_dir}/results", mode: 'copy', pattern: "*.csv"
    publishDir "${params.out_dir}/run/analysis", mode: 'copy', pattern: "analysis.log"
    publishDir "${params.out_dir}/run/analysis", mode: 'copy', pattern: "best_designs.jsonl"

    input:
    path pdbs

    output:
    path "best_designs.jsonl", emit: jsonl, topic: metadata_ch_fold_seq
    path "analysis.log", emit: log
    path "relaxed_pdbs/*.pdb", emit: relaxed_pdbs

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
    publishDir "${params.out_dir}/run/filter_analysis", mode: 'copy', pattern: 'filtered.jsonl'

    input:
    path(jsonl_file)
    path(pdb_files)
    val(paramString)

    output:
    path ("output/*.pdb"), emit: pdbs, optional: true
    path "filter_analysis_${task.index}.log"
    path ("filtered.jsonl"), emit: jsonl, optional: true

    script:

    """    
    python -u /scripts/filter_analysis.py \
        --jsonl-file ${jsonl_file} \
        --pdb-directory ./ \
        ${paramString} \
        --output-directory output \
        2>&1 | tee filter_analysis_${task.index}.log
    """
}
