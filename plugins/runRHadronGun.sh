#!/bin/bash

# The root directory to start the iteration from
root_directory=$1

# Function to iterate through subdirectories and their files
iterate_directories() {
    local dir=$1
    # Print the name of the current directory
    echo "Directory: $dir"

    # Iterate through the files in the current directory
    for file in "$dir"/*; do
        if [ -f "$file" ]; then
            echo "File: $(basename "$file")"
        fi
    done

    # Recursively iterate through each subdirectory
    for sub_dir in "$dir"/*/; do
        if [ -d "$sub_dir" ]; then
            iterate_directories "$sub_dir"
        fi
    done
}


# Create a dictionary to store key-value pairs for R-Hadron name and PDGID
declare -A d
d["~g_rho+"]=1009213
my_dict["key2"]="value2"
my_dict["key3"]="value3"

# Check if the directory argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Start the iteration from the root directory
iterate_directories "$root_directory"
