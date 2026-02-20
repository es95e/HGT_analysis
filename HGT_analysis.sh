#!/bin/bash
#
CONDA_SH="/home/your_user/miniconda/etc/profile.d/conda.sh" # insert your user or your miniconda path

if [ ! -f "$CONDA_SH" ]; then
    echo " Error: Conda not found at $CONDA_SH"
    echo "Please check the path or modify the script."
    exit 1
fi

echo "Activating Conda environment 'cub'..."
source "$CONDA_SH"
conda activate cub

run_script() {
    local script="$1"
    echo "=============================="
    echo "Executing: $script"
    echo "=============================="
    python "$script"
    if [ $? -eq 0 ]; then
        echo " $script executed successfully."
    else
        echo " Error during the execution of $script."
        echo "Aborting pipeline."
        conda deactivate
        exit 1
    fi
    echo
}

run_script "cub_calculator.py"
run_script "HGT_identification.py"
run_script "HGT_extractor.py"

echo "Deactivating Conda environment 'cub'..."
conda deactivate
echo "Pipeline completed successfully."
