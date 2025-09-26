#!/bin/bash

# Quick Tests for BindSweeper - Tests noise, scaffolds, and hotspots sweeps
# Usage: ./quick_tests.sh [output_base_dir]

set -euo pipefail

# TEST CONFIGURATION - Set to true/false to enable/disable tests
ENABLE_NOISE_SWEEP=false    # Test 1: Noise Scale Sweep
ENABLE_SCAFFOLD_SWEEP=false  # Test 2: Scaffold Sweep  
ENABLE_HOTSPOTS_SWEEP=true # Test 3: Hotspots Sweep
ENABLE_MULTI_SWEEP=false # Test 4: Multi-dimensional Sweep (hotspots + noise)
ENABLE_FILTERING_TEST=false # Test 5: Filtering Parameters Test

# Show help if requested
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    echo "Usage: $0 [output_base_dir]"
    echo ""
    echo "Quick Tests for BindSweeper - Tests noise, scaffolds, and hotspots sweeps"
    echo ""
    echo "Arguments:"
    echo "  output_base_dir    Base directory for test outputs (default: ./quick_test_results)"
    echo ""
    echo "This script runs five quick tests to validate BindSweeper functionality:"
    echo "  1. Noise Scale Sweep - Tests different noise levels for binder generation"
    echo "  2. Scaffold Sweep - Tests different fold conditioning scaffolds"
    echo "  3. Hotspots Sweep - Tests different hotspot combinations"
    echo "  4. Multi-dimensional Sweep - Tests combinations of hotspots and noise levels"
    echo "  5. Filtering Parameters Test - Tests filtering parameter sweeps"
    echo ""
    echo "Test Configuration:"
    echo "  You can enable/disable individual tests by editing the configuration"
    echo "  variables at the top of this script:"
    echo "    ENABLE_NOISE_SWEEP=true/false"
    echo "    ENABLE_SCAFFOLD_SWEEP=true/false"
    echo "    ENABLE_HOTSPOTS_SWEEP=true/false"
    echo "    ENABLE_MULTI_SWEEP=true/false"
    echo "    ENABLE_FILTERING_TEST=true/false"
    echo ""
    echo "Examples:"
    echo "  $0                    # Use default output directory"
    echo "  $0 /tmp/my_tests      # Use custom output directory"
    exit 0
fi

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Check if we're in the right directory and required files exist
if [[ ! -f "test_configs/test_noise_sweep.yaml" ]] || [[ ! -f "test_configs/test_scaffold_sweep.yaml" ]] || [[ ! -f "test_configs/test_hotspots_sweep.yaml" ]] || [[ ! -f "test_configs/test_multi_sweep.yaml" ]] || [[ ! -f "test_configs/test_filtering.yaml" ]]; then
    echo -e "${RED}Error: Required test configuration files not found in test_configs directory${NC}"
    echo "Expected files: test_configs/test_noise_sweep.yaml, test_configs/test_scaffold_sweep.yaml, test_configs/test_hotspots_sweep.yaml, test_configs/test_multi_sweep.yaml, test_configs/test_filtering.yaml"
    echo "Current directory: $SCRIPT_DIR"
    exit 1
fi

if [[ ! -f "../nextflow.config" ]]; then
    echo -e "${RED}Error: nextflow.config not found in parent directory${NC}"
    echo "Expected location: $(realpath ../nextflow.config)"
    exit 1
fi

# Check if uv is available
if ! command -v uv &> /dev/null; then
    echo -e "${RED}Error: 'uv' command not found${NC}"
    echo "Please ensure uv is installed and available in your PATH"
    echo "See: https://docs.astral.sh/uv/getting-started/installation/"
    exit 1
fi

# Default output directory
OUTPUT_BASE_DIR="${1:-/vast/scratch/users/$USER/quicktest_results}"
OUTPUT_DIR="${OUTPUT_BASE_DIR}"

echo -e "${BLUE}=== BindSweeper Quick Tests ===${NC}"
echo "Output directory: ${OUTPUT_DIR}"
echo

# Function to run a test
run_test() {
    local test_name="$1"
    local config_file="$2"
    local test_output_dir="${OUTPUT_DIR}/${test_name}_quicktest"
    
    echo -e "${YELLOW}--- Running ${test_name} test ---${NC}"
    echo "Config: ${config_file}"
    echo "Output: ${test_output_dir}"
    
    # Run the test with quick test mode (reduces designs to 2 each)
    if uv run bindsweeper --config "${config_file}" --output-dir "${test_output_dir}" --nextflow-config "../nextflow.config" --pipeline-path ../main.nf --quick-test --dry-run; then
        echo -e "${GREEN}✓ ${test_name} test configuration validated successfully${NC}"
        
        # Run the test automatically without user input
        echo -e "${YELLOW}Running ${test_name} test (this may take a while)...${NC}"
        if uv run bindsweeper --config "${config_file}" --output-dir "${test_output_dir}" --nextflow-config "../nextflow.config" --pipeline-path ../main.nf --quick-test --yes-to-all; then
            echo -e "${GREEN}✓ ${test_name} test completed successfully${NC}"
        else
            echo -e "${RED}✗ ${test_name} test failed${NC}"
            return 1
        fi
    else
        echo -e "${RED}✗ ${test_name} test configuration failed${NC}"
        return 1
    fi
    echo
}

# Ensure output directory exists
mkdir -p "${OUTPUT_DIR}"

echo -e "${BLUE}Running 5 quick tests to validate BindSweeper functionality:${NC}"
echo "1. Noise Scale Sweep - Tests different noise levels for binder generation"
echo "2. Scaffold Sweep - Tests different fold conditioning scaffolds" 
echo "3. Hotspots Sweep - Tests different hotspot combinations"
echo "4. Multi-dimensional Sweep - Tests combinations of hotspots and noise levels"
echo "5. Filtering Parameters Test - Tests filtering parameter sweeps"
echo

# Test 1: Noise Scale Sweep
if [[ "$ENABLE_NOISE_SWEEP" == "true" ]]; then
    run_test "noise_sweep" "test_configs/test_noise_sweep.yaml"
else
    echo -e "${YELLOW}⚠ Skipping noise sweep test (disabled)${NC}"
fi

# Test 2: Scaffold Sweep  
if [[ "$ENABLE_SCAFFOLD_SWEEP" == "true" ]]; then
    run_test "scaffold_sweep" "test_configs/test_scaffold_sweep.yaml"
else
    echo -e "${YELLOW}⚠ Skipping scaffold sweep test (disabled)${NC}"
fi

# Test 3: Hotspots Sweep
if [[ "$ENABLE_HOTSPOTS_SWEEP" == "true" ]]; then
    run_test "hotspots_sweep" "test_configs/test_hotspots_sweep.yaml"
else
    echo -e "${YELLOW}⚠ Skipping hotspots sweep test (disabled)${NC}"
fi

# Test 4: Multi-dimensional Sweep
if [[ "$ENABLE_MULTI_SWEEP" == "true" ]]; then
    run_test "multi_sweep" "test_configs/test_multi_sweep.yaml"
else
    echo -e "${YELLOW}⚠ Skipping multi-dimensional sweep test (disabled)${NC}"
fi

# Test 5: Filtering Parameters Test
if [[ "$ENABLE_FILTERING_TEST" == "true" ]]; then
    run_test "filtering_test" "test_configs/test_filtering.yaml"
else
    echo -e "${YELLOW}⚠ Skipping filtering parameters test (disabled)${NC}"
fi

echo -e "${GREEN}=== All quick tests completed! ===${NC}"
echo "Results saved to: ${OUTPUT_DIR}"
echo
echo -e "${BLUE}Summary of tests:${NC}"
echo "• Noise sweep: Tests rfd_noise_scale parameter (0.0, 0.2, 0.4)"
echo "• Scaffold sweep: Tests 3 different scaffold directories"  
echo "• Hotspots sweep: Tests different hotspot combinations"
echo "• Multi-dimensional sweep: Tests 2 hotspot combinations × 2 noise levels (4 total combinations)"
echo "• Filtering parameters test: Tests different filtering parameter combinations"
