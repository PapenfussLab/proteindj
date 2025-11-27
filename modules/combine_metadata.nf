process CombineMetadata {
    label 'python_tools'
    cache false
    publishDir "${params.out_dir}/run/combine_metadata", mode: 'copy'

    input:
    path metadata_fold
    path metadata_fold_seq

    output:
    path ("combined_metadata.csv"), emit: csv
    path "combined_metadata.log"

    script:
    """
    #!/bin/bash
    
    python /scripts/metadata_converter.py \
        --input_file ${metadata_fold} ${metadata_fold_seq} \
        --type csv \
        --output_file combined_metadata.csv \
        2>&1 | tee combined_metadata.log
    """
}
