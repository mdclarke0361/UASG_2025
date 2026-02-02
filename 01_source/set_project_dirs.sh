#!/bin/bash

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
