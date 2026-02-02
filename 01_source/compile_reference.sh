#!/bin/bash

# Compile VCF Reference File

# Source initialization script
source 01_source/initialize_script.sh

# Read in arguments
genotypes_normalized_file=${project_dir}/${1}
ref_vcf_dir=${project_dir}/${2}

# Output files
output_file=${reference_data_dir}/concatenated_ref.vcf.gz

# Create file list for concat - must be sorted correctly by chromosome
ref_vcf_file_list=${reference_data_dir}/ref_vcf_file_list.txt

for file in "${ref_vcf_dir}"/*.gz; do 
	echo $file
done | sort -V > $ref_vcf_file_list

# Run concatenation
bcftools concat \
	--output $output_file \
	--output-type z \
	--threads $threads \
	--file-list $ref_vcf_file_list \
	> $log_file

# Index the new concatenated ref file.
bcftools index \
	--threads $threads \
	$output_file

# test the chromosome names of our bead chip data against the chromosome names in the reference. 
# They must match for the imputation software to run.
bcftools query \
	-f '%CHROM\n' \
	$genotypes_normalized_file \
	| uniq \
	> ${report_out}/genotypes_chr_list.txt

bcftools query \
	-f '%CHROM\n' \
	$output_file \
	| uniq \
	> ${report_out}/ref_vcf_chr_list.txt

# Report any differences in chromosome numbers
diff -y \
	${report_out}/genotypes_chr_list.txt \
	${report_out}/ref_vcf_chr_list.txt \
	| tee -a $log_file
