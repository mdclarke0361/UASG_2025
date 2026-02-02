#!/bin/bash

# Phasing and Imputation of Missing Genotypes with BEAGLE 5.5

# Source initialization script
source 01_source/initialize_script.sh

# Read in arguments
memory_to_use=${1}
genotypes_fixed_ploidy_file=${project_dir}/${2}
concatenated_ref_vcf=${project_dir}/${3}
concatenated_map_file=${project_dir}/${4}
dr2_threshold=${5}

# Set output file names
genotypes_imputed=${processed_data_dir}/genotypes_imputed.vcf.gz
genotypes_imputed_filtered=${processed_data_dir}/genotypes_imputed_filtered.bcf

# Strip extension from output filename for use with BEAGLE
genotypes_imputed_basename=${genotypes_imputed%%.*}

# Run imputation software
beagle \
	-Xmx${memory_to_use}g \
	gt=${genotypes_fixed_ploidy_file} \
	ref=${concatenated_ref_vcf} \
	out=${genotypes_imputed_basename} \
	map=${concatenated_map_file} \
	nthreads=$threads

# Send the software created log file to log dir
mv ${genotypes_imputed_basename}.log $log_file

# Index the imputed file
bcftools index \
	--threads $threads \
	$genotypes_imputed

# Filter the imputed data

# Create expression to direct bcftools to filter on given threshold
# filtering_expression='INFO/DR2<"${dr2_threshold}"'

bcftools view \
	--exclude "INFO/DR2<${dr2_threshold}" \
	-O u \
	-o $genotypes_imputed_filtered \
	--threads $threads \
	$genotypes_imputed

# Index the filtered files
bcftools index \
	--threads $threads \
	$genotypes_imputed_filtered

# Report to stdout the output file name and path.
echo -e "${GREEN}Imputation complete... Output file is located at : ${YELLOW}$genotypes_imputed_filtered${NC}"
