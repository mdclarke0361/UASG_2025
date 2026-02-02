#!/bin/bash

# SNP Processing Step 1. Genotype Calling

# Source initialization script
source 01_source/initialize_script.sh

# Read in arguments
bpm_manifest=${project_dir}/${1}
cluster_file=${project_dir}/${2}
idat_dir=${project_dir}/${3}
sample_manifest=${project_dir}/${4}
csv_manifest=${project_dir}/${5}
ref_genome_fasta=${project_dir}/${6}

# Set output filenames and paths
gtc_out_dir=${processed_data_dir}/gtc_files
mkdir -p $gtc_out_dir

vcf_out_dir=${processed_data_dir}/vcf_files
mkdir -p $vcf_out_dir

output_file=${processed_data_dir}/genotypes.bcf

# Add the relative path for the DRAGEN binary to PATH variable. 
# This will allow execution of the 'dragena' command.
export PATH="${project_dir}/01_source/dragena:$PATH"

# Convert intensities to GTC files
echo -e "${YELLOW}Converting intensities to genotype calls in gtc format...${NC}"

# Run DRAGEN
dragena genotype call \
	--bpm-manifest $bpm_manifest \
	--cluster-file $cluster_file \
	--idat-folder $idat_dir \
	--sample-sheet $sample_manifest \
	--num-threads $threads \
	--output-folder $gtc_out_dir \
	> $log_file

# Run the file conversion using DRAGEN
echo -e "${YELLOW}Converting gtc files to vcf files...${NC}"

# Run DRAGEN
dragena genotype gtc-to-vcf \
	--bpm-manifest $bpm_manifest \
	--csv-manifest $csv_manifest \
	--genome-fasta-file $ref_genome_fasta \
	--gtc-folder $gtc_out_dir \
	--sample-sheet $sample_manifest \
	--output-folder $vcf_out_dir \
	>> $log_file

# Merge and index sample files
echo -e "${YELLOW}Merging vcf files and outputting in bcf format...${NC}"

# Create file list for VCF files
vcf_file_arr=()

for file in "${vcf_out_dir}"/*.snv.vcf.gz; do

	vcf_file_arr+=("$file")

done

vcf_file_list="${vcf_file_arr[*]}"

# Index the files
for file in $vcf_file_list; do

	bcftools index \
		--threads $threads \
		$file

done

# Merge the files and output as a BCF. Don't include a lengthy addition to the header.
bcftools merge \
	--output $output_file \
	--output-type u \
	--threads $threads \
	--no-version \
	$vcf_file_list \
	>> $log_file

# Report to stdout the output file name and path.
echo -e "${GREEN}Genotype calling complete... Output file is located at : ${YELLOW}$output_file${NC}"