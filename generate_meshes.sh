#!/bin/bash

# List of resolutions (edit as needed)
resolutions=(10 20 30 40 50)

# Base command path
BIN_DIR="./bin"
OUTPUT_DIR="./output"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Mesh names
meshes=("COMesh" "CSMesh")

for mesh in "${meshes[@]}"; do
    for res in "${resolutions[@]}"; do
        output_file="${OUTPUT_DIR}/${mesh}_${res}.g"
        
        if [ "$mesh" == "COMesh" ]; then
            echo "Generating ICOMesh at resolution $res -> $output_file"
            "$BIN_DIR/GenerateICOMesh" --res "$res" --dual --file "$output_file"
        elif [ "$mesh" == "CSMesh" ]; then
            echo "Generating CSMesh at resolution $res -> $output_file"
            "$BIN_DIR/GenerateCSMesh" --res "$res" --alt --file "$output_file"
        else
            echo "Unknown mesh type: $mesh"
        fi
    done
done
