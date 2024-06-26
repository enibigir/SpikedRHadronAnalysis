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
d["~g_K*0"]=1009313
d["~g_K*+"]=1009323
d["~g_rho0"]=1009113
d["~g_omega"]=1009223
d["~g_phi"]=1009333
d["~g_Delta-"]=1091114
d["~g_Delta0"]=1092114
d["~g_Delta+"]=1092214
d["~g_Delta++"]=1092224
d["~g_Sigma*-"]=1093114
d["~g_Sigma*0"]=1093214
d["~g_Sigma*+"]=1093224
d["~g_Xi*-"]=1093314
d["~g_Xi*0"]=1093324
d["~g_Omega-"]=1093334
d["~g_rho-"]=-1009213
d["~g_K*bar0"]=-1009313
d["~g_K*-"]=-1009323
d["~g_Deltabar+"]=-1091114
d["~g_Deltabar0"]=-1092114
d["~g_Deltabar-"]=-1092214
d["~g_Deltabar--"]=-1092224
d["~g_Sigma*bar+"]=-1093114
d["~g_Sigma*bar0"]=-1093214
d["~g_Sigma*bar-"]=-1093224
d["~g_Xi*bar+"]=-1093314
d["~g_Xi*bar0"]=-1093324
d["~g_Omegabar+"]=-1093334

# Check if the directory argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Start the iteration from the root directory
iterate_directories "$root_directory"
