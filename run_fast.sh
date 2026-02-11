#!/bin/bash
# Check if sysimage exists
if [ -f "pmd_sysimage.so" ]; then
    julia --project=. --sysimage pmd_sysimage.so "$@"
else
    echo "Warning: Sysimage 'pmd_sysimage.so' not found. Running slow..."
    julia --project=. "$@"
fi
