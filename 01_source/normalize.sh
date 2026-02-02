#!/bin/bash

# Normalization of Genotype Calls

# Source initialization script
source 01_source/initialize_script.sh

# Read in arguments
genotypes_filtered_controls_file=${project_dir}/${1}
ref_genome_fasta=${project_dir}/${2}

# Output files
output_file=${processed_data_dir}/genotypes_normalized.bcf

# Normalize Genotype Calls with bcftools
bcftools norm \
	--check-ref w \
	--fasta-ref $ref_genome_fasta \
	--output-type u \
	--output $output_file \
	--threads $threads \
	$genotypes_filtered_controls_file \
	2> $log_file

# Index the normalized BCF
bcftools index \
	--threads $threads \
	$output_file

# Print the resulting output file path
echo -e "${GREEN}Normalization complete...  Output file is located at : ${YELLOW}$output_file${NC}"
