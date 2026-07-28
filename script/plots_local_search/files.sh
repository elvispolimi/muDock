#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_folder> <output_file>"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_FILE="$2"

# Empty the output file if it already exists
> "$OUTPUT_FILE"

for file in "$INPUT_DIR"/*; do
    if [[ -f "$file" ]]; then
        cat "$file" >> "$OUTPUT_FILE"
        echo >> "$OUTPUT_FILE"   # Add a newline between files
    fi
done

echo "Concatenated files from '$INPUT_DIR' into '$OUTPUT_FILE'."