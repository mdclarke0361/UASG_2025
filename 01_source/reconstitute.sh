#!/bin/bash

# <PROGRAM TITLE>

# Source initialization script
source 01_source/initialize_script.sh

# Read in arguments
genotypes_imputed_filtered_file=${project_dir}/${1}
genotypes_normalized_file=${project_dir}/${2}

# Output files
genotypes_imputed_only=${processed_data_dir}/genotypes_imputed_only.bcf
genotypes_combined=${processed_data_dir}/genotypes_combined.bcf

output_file=${processed_data_dir}/genotypes_reconstituted.bcf

# Get the number of threads to use based on desired processor use.
core_total="$(nproc --all)"
threads=$(echo "${core_total}*0.8" | bc -l)
threads=${threads%.*}

# Separate out the imputed markers from the imputed data
bcftools view \
	--include 'INFO/IMP=1' \
	-O u \
	-o $genotypes_imputed_only \
	--threads $threads \
	$genotypes_imputed_filtered_file

bcftools index \
	--threads $threads \
	$genotypes_imputed_only

# Test that all measured markers have been removed.
measured_marker_count=$(bcftools view -H $genotypes_imputed_only | grep -cv ";IMP")

echo -e "${YELLOW}$measured_marker_count${NC} measured markers remaining in the imputed data..." \
	| tee $log_file

bcftools concat \
	--allow-overlaps \
	--output $genotypes_combined \
	--output-type u \
	--threads $threads \
	$genotypes_imputed_only \
	$genotypes_normalized_file

# Sort the files to return to positional order.
bcftools sort \
	--output $output_file \
	-O u \
	-v 4 \
	$genotypes_combined

# Index the sorted file.
bcftools index \
	--min-shift 30 \
	--threads $threads \
	$output_file

# Report to stdout the output file name and path.
echo -e "${GREEN}Genotype calling complete... Output file is located at : ${YELLOW}$output_file${NC}"
