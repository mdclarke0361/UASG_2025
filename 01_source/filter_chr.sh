#!/bin/bash

# Filter out Chromosomes Unable to be Imputed

# Source initialization script
source 01_source/initialize_script.sh

# Read in arguments
genotypes_normalized_file=${project_dir}/${1}
chr_list_file=${project_dir}/${2}
chr_to_remove=${3}

output_file=${processed_data_dir}/genotypes_filtered_chr.vcf.gz

# Get the list of chromosomes in the genotype call file
chr_list=$(<$chr_list_file)

chr_to_remove=${chr_to_remove/','/' '}

# Remove chromosomes based on a given list.
for chr in $chr_to_remove; do

	chr_list=${chr_list/$chr/}

done

# Transform into comma-separated list.
chr_list=$(echo $chr_list | tr " " ",")

echo -e "${YELLOW}The following chromosomes will be kept: $chr_list${NC}" \
	| tee $log_file

# Run the filtering
bcftools view \
	--output $output_file \
	--output-type z \
	--regions $chr_list \
	--threads $threads \
	$genotypes_normalized_file

# Test the output
kept_chr_list=$(
	bcftools query \
	-f '%CHROM\n' \
	$output_file \
	| uniq
)

echo -e "${GREEN}Incompatible chromosomes removed... \
	The following chromosomes remain: ${YELLOW}$kept_chr_list${NC}" \
	| tee -a $log_file

echo -e "${GREEN}Output file located at $output_file"