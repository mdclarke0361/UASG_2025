#!/bin/bash

# Filter Out Controls

# Source initialization script
source 01_source/initialize_script.sh

# Read in arguments
genotypes_fixed_header_file=${project_dir}/${1}
sample_list_file=${project_dir}/${2}
control_sample_prefix=${3}

#
output_file=${processed_data_dir}/genotypes_filtered_controls.bcf

# Create a comma-separated exclusion list for control samples.
exclude_list=$(grep "^${control_sample_prefix}" $sample_list_file)
exclude_list="${exclude_list//$'\n'/','}"

#
bcftools view \
	--output $output_file \
	--samples ^$exclude_list \
	$genotypes_fixed_header_file \
	> $log_file

# Check the new list of samples
remaining_samples=$(
	bcftools query \
	--list-samples \
	$output_file \
	| wc -l
)

# Print the number of remaining (non-control) samples
echo -e "${YELLOW}There are ${remaining_samples} remaining samples${NC}" | tee -a $log_file

# Print file path of output file
echo -e "${GREEN}Filtering complete...  Output file is located at : ${YELLOW}$output_file${NC}"