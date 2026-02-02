#!/bin/bash

# Rename the Samples in File Header

# Source initialization script
source 01_source/initialize_script.sh

# Read in arguments
genotypes_file=${project_dir}/${1}
sample_manifest=${project_dir}/${2}

#
output_file=${processed_data_dir}/genotypes_fixed_header.bcf

# Get the list of sample names used in the file header and save them to a file.
barcode_list_file=${metadata_dir}/snp_sample_barcodes.txt

bcftools query \
	--list-samples \
	$genotypes_file \
	> $barcode_list_file

# Read the barcodes line by line and get the corresponding sample name from the sample list.
sample_list_file=${metadata_dir}/snp_sample_names.txt
touch $sample_list_file

while read -r line; do

	barcode=${line%_*}
	position=${line#*_}

	sample_name=$(
		grep $barcode $sample_manifest | 
		grep $position | 
		cut -d "," -f 1
	)

	# replace the hyphen with an underscore
	sample_name=${sample_name/-/_}

	echo $sample_name >> $sample_list_file

done < $barcode_list_file

#
bcftools reheader \
	--samples $sample_list_file \
	--output $output_file \
	--threads $threads \
	$genotypes_file \
	> $log_file

# Report to stdout the output file name and path.
echo -e "${GREEN}Header correction complete... Output file is located at : ${YELLOW}$output_file${NC}"