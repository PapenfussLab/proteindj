process Compress {
    label 'pyrosetta_tools'
    publishDir "${params.out_dir}/run/${program}", mode: 'copy', pattern: "*.tar.gz"

    input:
    val program
    // program which produced PDB files e.g. rfd, af2, boltz, fampnn, mpnn
    path files

    output:
    path "*.tar.gz"

    script:
    def num_processes = task.cpus - 1
    """
    # To avoid line length limits, first need to copy files into a new directory
    # Will exclude nextflow files .command.out, .command.err etc.
    mkdir ${program}_results/
    for file in \$(find *.*); do
        cp -L "\$file" ./${program}_results/
    done

    # Uses pigz for parallel compression
    tar -h \
        --use-compress-program="pigz -p ${num_processes}" \
        -cf ${program}_results.tar.gz \
        ${program}_results
    """
}
