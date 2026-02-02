#!/bin/bash

# Convert Haploid X Chromosome Markers to Diploid

# Source initialization script
source 01_source/initialize_script.sh

# Read in arguments
genotypes_filtered_chr_file=${project_dir}/${1}
chr_length=${2}

output_file=${processed_data_dir}/genotypes_fixed_ploidy.vcf.gz

#
ploidy_temp_file=${reference_data_dir}/ploidy.txt

echo "X 1 $chr_length M 2" > $ploidy_temp_file

#
bcftools +fixploidy \
	--output-type z \
	--output $output_file \
	--threads $threads \
	$genotypes_filtered_chr_file \
	-- \
	--ploidy $ploidy_temp_file \
	> $log_file

rm $ploidy_temp_file

echo -e "${GREEN}X chromosome corrected... output file located at: $output_file${NC}"