process AnalyseBestDesigns {
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
