#!/bin/bash

# SEATREE Python3 Cache Cleanup Script
# Removes Python cache files and bytecode

set -e

echo "SEATREE Python3 Cache Cleanup"
echo "=============================="

# Function to safely remove files/directories
safe_remove() {
    local pattern="$1"
    local description="$2"
    local count=$(find . -name "$pattern" | wc -l)

    if [ "$count" -gt 0 ]; then
        echo "Removing $count $description..."
        find . -name "$pattern" -delete
    else
        echo "No $description found"
    fi
}

echo ""
echo "Cleaning Python cache files..."

# Remove Python bytecode files first (so directories become empty)
safe_remove "*.pyc" ".pyc files"
safe_remove "*.pyo" ".pyo files"
safe_remove "*.pyd" ".pyd files"

# Now remove empty __pycache__ directories
echo "Removing __pycache__ directories..."
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true

# Remove other Python temporary files
safe_remove ".coverage" "coverage files"

echo ""
echo "Python cache cleanup complete!"
echo ""
echo "What was cleaned:"
echo "  ✓ __pycache__ directories"
echo "  ✓ .pyc/.pyo bytecode files"
echo "  ✓ Other Python temporary files"