#!/bin/bash

# upload_containers.sh
# Usage: ./upload_containers.sh v2.2

set -e  # Exit on error

# Configuration
REGISTRY="oras://ghcr.io/papenfusslab/proteindj"
CONTAINERS=(
    "af2"
    "bindcraft"
    "boltz2"
    "dl_binder_design"
    "fampnn"
    "pyrosetta_tools"
    "python_tools"
    "rfdiffusion"
)

# Check if version argument provided
if [ -z "$1" ]; then
    echo "Error: No version specified"
    echo "Usage: $0 <version>"
    echo "Example: $0 v2.2"
    exit 1
fi

VERSION="$1"

# Validate version format (optional)
if [[ ! $VERSION =~ ^v[0-9]+\.[0-9]+$ ]]; then
    echo "Warning: Version format doesn't match vX.Y pattern (e.g., v2.2)"
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

echo "========================================="
echo "Uploading containers with tags:"
echo "  - Version: $VERSION"
echo "  - latest"
echo "========================================="
echo ""

# Counter for progress
TOTAL=${#CONTAINERS[@]}
COUNT=0

# Loop through each container
for CONTAINER in "${CONTAINERS[@]}"; do
    COUNT=$((COUNT + 1))
    SIF_FILE="${CONTAINER}.sif"
    
    echo "[$COUNT/$TOTAL] Processing $SIF_FILE..."
    
    # Check if SIF file exists
    if [ ! -f "$SIF_FILE" ]; then
        echo "  ERROR: $SIF_FILE not found - skipping"
        continue
    fi
    
    # Push with version tag
    echo "  Pushing ${CONTAINER}:${VERSION}..."
    if apptainer push "$SIF_FILE" "${REGISTRY}/${CONTAINER}:${VERSION}"; then
        echo "  ✓ Successfully pushed ${CONTAINER}:${VERSION}"
    else
        echo "  ✗ Failed to push ${CONTAINER}:${VERSION}"
        continue
    fi
    
    # Push with latest tag
    echo "  Pushing ${CONTAINER}:latest..."
    if apptainer push "$SIF_FILE" "${REGISTRY}/${CONTAINER}:latest"; then
        echo "  ✓ Successfully pushed ${CONTAINER}:latest"
    else
        echo "  ✗ Failed to push ${CONTAINER}:latest"
    fi
    
    echo ""
done

echo "========================================="
echo "Upload complete!"
echo "========================================="
echo ""
echo "Containers available at:"
for CONTAINER in "${CONTAINERS[@]}"; do
    echo "  - ${REGISTRY}/${CONTAINER}:${VERSION}"
    echo "  - ${REGISTRY}/${CONTAINER}:latest"
done
