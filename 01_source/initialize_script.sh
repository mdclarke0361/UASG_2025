#!/bin/bash

# Initialize project scripts.

# Ensure all project directories exist, if not, create them.

subdirectories=(
	# Data directories
	"02_data/metadata" "02_data/processed" "02_data/raw" "02_data/reference"
	# Results directories
	"03_results/figures" "03_results/reports" "03_results/tables" "03_results/logs"
)

for dir in "${subdirectories[@]}"; do

	mkdir -p $dir

done

# Set up of text colors for terminal output.
# Text Colors
export RED='\033[0;31m'
export GREEN='\033[0;32m'
export YELLOW='\033[0;33m'
export NC='\033[0m' # Reset to default

# Get path for project directory
export metadata_dir="02_data/metadata"
export processed_data_dir="02_data/processed"
export raw_data_dir="02_data/raw"
export reference_data_dir="02_data/reference"
export figure_out="03_results/figures"
export report_out="03_results/reports"
export table_out="03_results/tables"
export log_out="03_results/logs"

# Get the number of threads to use based on desired processor use.
core_total="$(nproc --all)"
threads=$(echo "${core_total}*0.8" | bc -l)
export threads=${threads%.*}

# Get path for project directory
script_name=$(basename $0)
script_dir=$(dirname $0)
export project_dir=$(dirname $script_dir)

# Assign log filename
export log_file="${log_out}/${script_name%%.*}.log"
