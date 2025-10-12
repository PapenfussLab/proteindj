#!/bin/bash

# ProteinDJ Pipeline Testing Framework
# ===========================================
# This script tests all pipeline modes and generates a comprehensive report

set -uo pipefail  # Exit on undefined vars and pipe failures (but allow command failures)

# Configuration
COMPUTE="$1"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

OUT_DIR="$2"/"${TIMESTAMP}"
LOG_DIR="${OUT_DIR}/logs"

# Create directories
mkdir -p "${LOG_DIR}"

# Test modes to run
declare -a MODES=(
    "monomer_denovo"
    "monomer_foldcond" 
    "monomer_motifscaff"
    "monomer_partialdiff"
    "binder_denovo"
    "binder_foldcond"
    "binder_motifscaff"
    "binder_partialdiff"
)

# Function to log with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_DIR}/test_execution.log"
}

# Function to run a single test mode
run_test_mode() {
    local mode=$1
    local start_time=$(date +%s)
    local log_file="${LOG_DIR}/${mode}.log"
    local status="UNKNOWN"
    local error_msg=""
    local output_dir=""
    
    log "Starting test for mode: ${mode}, compute profile: ${COMPUTE}"
    
    # Set mode-specific output directory
    output_dir="${OUT_DIR}/${mode}"
    
    # Run the pipeline
    if nextflow run main.nf \
        -profile "${COMPUTE},test,${mode}" \
        --zip_pdbs false \
        --out_dir "${output_dir}" \
        > "${log_file}" 2>&1; then
        status="PASSED"
        log "✅ ${mode} test PASSED"
    else
        status="FAILED"
        error_msg=$(tail -n 20 "${log_file}" | grep -E "(ERROR|Exception|Failed)" | head -n 3 || echo "Check log file for details")
        log "❌ ${mode} test FAILED"
    fi
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    # Collect test results
    echo "{
        \"mode\": \"${mode}\",
        \"status\": \"${status}\",
        \"duration_seconds\": ${duration},
        \"start_time\": \"$(date -d @${start_time} '+%Y-%m-%d %H:%M:%S')\",
        \"end_time\": \"$(date -d @${end_time} '+%Y-%m-%d %H:%M:%S')\",
        \"log_file\": \"${log_file}\",
        \"output_dir\": \"${output_dir}\",
        \"error_message\": \"${error_msg}\"
    }" >> "${OUT_DIR}/temp_results.json"
    
    return $([ "$status" = "PASSED" ] && echo 0 || echo 1)
}

# Function to analyze pipeline outputs
analyze_outputs() {
    local mode=$1
    local analysis_file="${OUT_DIR}/analysis.txt"
    # Set mode-specific output directory
    output_dir="${OUT_DIR}/${mode}"

    if [ -d "${output_dir}" ]; then
        echo "Output Analysis for ${mode}:" >> "${analysis_file}"
        echo "=================================" >> "${analysis_file}"
        echo "Intermediate output directories generated:" >> "${analysis_file}"
        basename -a ${output_dir}/run/*/ >> "${analysis_file}"
        find "${output_dir}" -type f -name "*.tar.gz*" | wc -l | xargs -I {} echo "Intermediate tar.gz archives generated: {}" >> "${analysis_file}"
        
        # Check for presence of both CSV files
        if [[ -f "${output_dir}/results/all_designs.csv" && -f "${output_dir}/results/best_designs.csv" ]]; then
            echo "Both all_designs.csv and best_designs.csv exist ✅" >> "${analysis_file}"
            
            lines_all=$(wc -l < "${output_dir}/results/all_designs.csv")
            echo "Number of final designs (all_designs.csv): $((lines_all - 1))" >> "${analysis_file}"
            
            lines_best=$(wc -l < "${output_dir}/results/best_designs.csv")
            echo "Number of best designs (best_designs.csv): $((lines_best - 1))" >> "${analysis_file}"
            
            pdb_count=$(find "${output_dir}/results/best_designs" -type f -name "*.pdb" | wc -l)
            echo "PDB files generated: ${pdb_count}" >> "${analysis_file}"
        else
            echo "Missing one or more CSV files (all_designs.csv, best_designs.csv) ❌" >> "${analysis_file}"
        fi
        
        # Directory size
        size=$(du -sh "${output_dir}" 2>/dev/null | cut -f1)
        echo "Total output size: ${size}" >> "${analysis_file}"
        echo "=================================" >> "${analysis_file}"
    else
        echo "Output directory not found: ${output_dir}" >> "${analysis_file}"
    fi
}

# Main execution
main() {
    log "🚀 Starting ProteinDJ Pipeline Test Suite"
    log "Test directory: ${OUT_DIR}"
    
    local failed_count=0
    local passed_count=0
    
    # Run tests for each mode
    for mode in "${MODES[@]}"; do
        if run_test_mode "${mode}"; then
            ((passed_count++))
            # Analyze outputs for successful runs
            analyze_outputs "${mode}"
        else
            ((failed_count++))
        fi
        
        # Add spacing between tests
        log "----------------------------------------"
    done
    
    # Final summary
    local analysis_file="${OUT_DIR}/analysis.txt"
    touch "${analysis_file}"
    log  "📊 Test Summary:"
    echo "📊 Test Summary:" >> "${analysis_file}"
    log  "   Passed: ${passed_count}/${#MODES[@]}"
    echo "   Passed: ${passed_count}/${#MODES[@]}" >> "${analysis_file}"
    log  "   Failed: ${failed_count}/${#MODES[@]}"
    echo "   Failed: ${failed_count}/${#MODES[@]}" >> "${analysis_file}"
    log  "   Success Rate: $((passed_count * 100 / ${#MODES[@]}))%"
    echo "   Success Rate: $((passed_count * 100 / ${#MODES[@]}))%" >> "${analysis_file}"
    log ""
    log  "📁 Results saved to: ${OUT_DIR}/"
    echo "📁 Results saved to: ${OUT_DIR}/" >> "${analysis_file}"
    log  "🔍 Logs: ${LOG_DIR}/"
    echo "🔍 Logs: ${LOG_DIR}/" >> "${analysis_file}"
    
    # Exit with appropriate code
    if [ ${failed_count} -eq 0 ]; then
        log "🎉 All tests passed!"
        exit 0
    else
        log "⚠️  Some tests failed. Check the report for details."
        exit 1
    fi
}

# Run main function
main "$@"

