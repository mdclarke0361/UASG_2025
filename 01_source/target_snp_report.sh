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
target_snp_list=${2}
processed_genotype_file=${project_dir}/${3}
ref_snp_list=${project_dir}/${4}

# Create directories for output
# Log files
log_file_dir=${analysis_dir}/log_files
mkdir -p $log_file_dir 2> /dev/null

log_file=${log_file_dir}/target_snp_report.log

# Output files
output_file=${analysis_dir}/target_snp_report.tsv

# Get the number of available CPU cores and calculate 80% of that number for threading.
core_total="$(nproc --all)"
threads=$(echo "${core_total}*0.8" | bc -l)
threads=${threads%.*}

# Create report, extract header 
bcftools view \
	-h $processed_genotype_file |
	tail -n 1 \
	> ${output_file}

# Begin counting any unrecognized snp IDs from target list
not_found=0

# Read the target SNP list line by line.
while read -r line; do

	# Strip any whitespace from the SNP name.
	target=${line//" "/""}

	echo -e "Searching for target SNP: ${YELLOW}$target${NC}" \
		| tee $log_file

	# Use grep to search for exact term and split line into elements of array.
	mapfile -t -d $'\t' snp_line \
		< <(grep -w $target $ref_snp_list) # grep -w for exact word match

	# Assign variables from array elements.
	chr=${snp_line[0]}
	start=${snp_line[1]}
	end=${snp_line[2]}
	ref_ID=${snp_line[3]//$'\n'/""}

	# Report the results of searching through reference SNP list.
	# If found, check genotype data for matching SNP ID at the position.
	if [[ -n $ref_ID ]]; then

		start_pos="$(( start - 1 ))"
		end_pos="$(( end + 1 ))"
		position="$chr:$start_pos-$end_pos"

		echo -e "\t${GREEN}Match found in reference database.${NC}" \
			| tee -a $log_file

		mapfile -t -d " " marker_arr \
			< <(bcftools query \
				--regions $position \
				--format '%CHROM %POS %ID' \
				$processed_genotype_file
			)

		rsID=${marker_arr[2]//$'\n'/""}

		# Compare reference SNP ID to genotype SNP ID at the position.
		if [[ -z $rsID ]]; then

			echo -e "\t${RED}Match not found in genotype data at position: $position.${NC}\n" \
				| tee -a  $log_file

			((not_found++))
			break

		elif [[ "$rsID" == "$ref_ID" ]]; then

			echo -e "\t${GREEN}Match found in genotype data.${NC}\n" \
				| tee -a $log_file

		else

			echo -e "\tMatch not found in genotype data." \
				"\tUsing ${YELLOW}$rsID ${NC}based on position.\n\tCorrecting rsID in report.\n" \
				| tee -a $log_file

		fi

		# Add line to report.
		# Split up any multiallelic markers (for ease of downstream analysis) and output to report.
		bcftools norm \
			--multiallelics - \
			--do-not-normalize \
			-r $position \
			$processed_genotype_file \
			2> $log_file |
			bcftools view \
			-H |
			sed "s/$rsID/$ref_ID/g" \
			>> ${output_file}
	
	else
	
		echo -e "\t${RED}Match not found in reference file.${NC}\n" \
			| tee -a $log_file
		((not_found++))

	fi

done < $target_snp_list

if [[ $not_found -gt 0 ]]; then

	echo -e "${RED}$not_found target SNPs cannot be found.${NC}\n" \
		| tee -a $log_file

	read -rp "Do you want to continue generating the report with the found genes? (y/n): " -n 1 -r
	

	if [[ $REPLY == "y" ]]; then

		echo -e "\n\n${GREEN}Continuing with report generation...${NC}\n" \
			| tee -a $log_file
	
	else

		echo -e "\n${RED}Report generation aborted by user.${NC}\n" \
			| tee -a $log_file

		exit 1

	fi

fi

# Report to stdout the output file name and path.
echo -e "${GREEN}Target SNP report complete... Output file is located at : ${YELLOW}$output_file${NC}" \
	| tee -a $log_file

# Call R script to clean the data.
Rscript ${script_dir}/clean_report.r \
	$output_file