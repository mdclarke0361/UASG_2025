#!/bin/bash

# Generate Genotype Report for Target SNPs

# Source initialization script
source 01_source/initialize_script.sh

# Read in arguments
target_snp_list=${1}
processed_genotype_file=${project_dir}/${2}
ref_snp_list=${project_dir}/${3}
report_file=${project_dir}/${4}

# Create report, extract header 
bcftools view \
	-h $processed_genotype_file |
	tail -n 1 \
	> ${report_file}

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
			>> ${report_file}
	
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
echo -e "${GREEN}Target SNP report complete... Output file is located at : ${YELLOW}$report_file${NC}" \
	| tee -a $log_file

# Call R script to clean the data.
Rscript ${script_dir}/clean_report.r \
	$report_file