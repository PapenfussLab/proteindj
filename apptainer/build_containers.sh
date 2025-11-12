#!/bin/bash

# ///  Configuration ///
BUILD_DIRECTORY="/dev/shm/${USER}/containers"

# Container build configuration
declare -A CONTAINERS=(
    ["af2"]="af2.def"
    ["bindcraft"]="bindcraft.def"
    ["boltz2"]="boltz2.def"
    ["dl_binder_design"]="dl_binder_design.def"
    ["fampnn"]="fampnn.def"
    ["pyrosetta_tools"]="pyrosetta_tools.def"
    ["rfdiffusion"]="rfdiffusion.def"
)

# Enable/disable container builds (1=enabled, 0=disabled)
BUILD_AF2=1
BUILD_BINDCRAFT=1
BUILD_BOLTZ2=1
BUILD_DL_BINDER_DESIGN=1
BUILD_FAMPNN=1
BUILD_PYROSETTA_TOOLS=1
BUILD_RFDIFFUSION=1

# Function to build a single container
build_container() {
    local name=$1
    local def_file=$2
    local job_name="build_${name}"
    local log_file="${BUILD_DIRECTORY}/${name}_build.log"
    
    # Create the build script for this container
    local build_script="${BUILD_DIRECTORY}/${name}_build.sh"
    
    cat << EOF > "$build_script"
#!/bin/bash
mkdir -p "${BUILD_DIRECTORY}"
export APPTAINER_TMPDIR="${BUILD_DIRECTORY}"

# Set up unbuffered output and redirect to log file
LOG_FILE="${log_file}"
exec 1> >(tee -a "\$LOG_FILE")
exec 2>&1

# NOTE: YOU MAY NEED TO ADJUST THIS FOR YOUR HPC ENVIRONMENT
module reset
module load apptainer
module load micromamba

echo "Starting build for ${name}..."
start_time=\$(date +%s)
apptainer build --fakeroot "${BUILD_DIRECTORY}/${name}.sif" ${def_file}
end_time=\$(date +%s)
duration=\$((end_time - start_time))
echo "${name} build completed in \$((duration / 60))m \$((duration % 60))s"
EOF

    chmod +x "$build_script"

    # Submit the build job
    # NOTE: YOU MAY NEED TO ADJUST THIS FOR YOUR HPC ENVIRONMENT
    sbatch \
        --job-name "$job_name" \
        --cpus-per-task 12 \
        --time 5:00:00 \
        --mem 48G \
        --output="/dev/null" \
        --error="/dev/null" \
        "$build_script"
        
    echo "Job submitted for ${name}. Log will be written to: ${log_file}"
}

# Create build directory if it doesn't exist
mkdir -p "${BUILD_DIRECTORY}"

# Submit build jobs for enabled containers
for container in "${!CONTAINERS[@]}"; do
    var_name="BUILD_${container^^}"
    var_name=${var_name//-/_}  # Replace hyphens with underscores
    
    if [[ "${!var_name}" == "1" ]]; then
        def_file="${CONTAINERS[$container]}"
        echo "Submitting build job for $container..."
        build_container "$container" "$def_file"
    else
        echo "Skipping disabled container: $container"
    fi
done

echo "All build jobs submitted. Check log files in ${BUILD_DIRECTORY} for progress."