#!/bin/bash

# Generate Genotype Report for Target SNPs

# Set up of text colors for terminal output.
# Text Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # Reset to default

# Get path for project directory
script_dir=$(dirname $0)
source_dir=$(dirname $script_dir)
project_dir=$(dirname $source_dir)

# Read in arguments
# Get required files from project subdirectories.
analysis_dir=${1}
target_gene_list=${2}
processed_genotype_file=${project_dir}/${3}
ref_gene_list=${project_dir}/${4}
bp_before_gene=${5}
bp_after_gene=${6}

# Create directories for output
# Log files
log_file_dir=${analysis_dir}/log_files
mkdir -p $log_file_dir 2> /dev/null

log_file=${log_file_dir}/target_gene_report.log

# Output files
output_file=${analysis_dir}/target_gene_report.tsv

# Get the number of available CPU cores and calculate 80% of that number for threading.
core_total="$(nproc --all)"
threads=$(echo "${core_total}*0.8" | bc -l)
threads=${threads%.*}

# Create report, extract header line from genotype file
bcftools view \
	-h $processed_genotype_file |
	tail -n 1 \
	> ${output_file}

# Instantiate array for gene locus information
gene_positions=()
not_found=0

# Read the target gene list line by line.
while read -r line; do

	# Strip any whitespace from the gene name.
	target=${line//" "/""}

	# Use grep to search for exact term and split line into elements of array.
	read -ra gene_line <<< "$(grep -w "$target" $ref_gene_list)"

	# Assign variables from array elements.
	chr=${gene_line[0]}
	start=${gene_line[1]}
	end=${gene_line[2]}
	gene=${gene_line[3]}

	# Report the results of searching through reference gene list.
	if [[ -n $gene ]]; then

		start_pos="$(( start - bp_before_gene ))"
		end_pos="$(( end + bp_after_gene ))"
		position="$chr:$start_pos-$end_pos"
		gene_positions+=("$position")

		echo -e "${GREEN}$gene found in reference.${NC}.\n" |
			tee $log_file

	else
	
		echo -e "${RED}$target not found in reference.${NC} Try another gene name.\n" |
			tee -a $log_file

		((not_found++))

	fi

done < $target_gene_list

if [[ $not_found -gt 0 ]]; then

	echo -e "${RED}$not_found target genes were not found in the reference gene list.${NC}\n" |
			tee -a $log_file

	read -rp "Do you want to continue generating the report with the found genes? (y/n): " -n 1 -r
	

	if [[ $REPLY == "y" ]]; then
		echo -e "\n${GREEN}Continuing with report generation...${NC}\n" |
			tee -a $log_file
	
	else
		echo -e "\n${RED}Report generation aborted by user.${NC}\n" |
			tee -a $log_file
			
		exit 1
	fi

fi

# Assemble target gene positions into a comma-separated list.
positions_list=$(echo ${gene_positions[*]} | tr " " ",")

# Split up any multiallelic markers (for ease of downstream analysis) and output to report.
bcftools norm \
	--multiallelics - \
	--do-not-normalize \
	--regions $positions_list \
	$processed_genotype_file \
	2> $log_file |
	bcftools view \
	-H \
	>> ${output_file}

# Report to stdout the output file name and path.
echo -e "${GREEN}Target gene report complete... Output file is located at : ${YELLOW}$output_file${NC}" \
	| tee -a $log_file

# Call R script to clean the data.
Rscript ${script_dir}/clean_report.r \
	$output_file