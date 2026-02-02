#!/bin/bash

# Gene Annotation for Genotype Data

# Source initialization script
source 01_source/initialize_script.sh

# Read in arguments
genotypes_reconstituted_file=${project_dir}/${1}
gene_list_file=${project_dir}/${2}

# 
output_file=${processed_data_dir}/genotypes_annotated.bcf

#
bcftools annotate \
	--annotations $gene_list_file \
	--columns CHROM,FROM,TO,GENE \
	--header-lines <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="GeneName">') \
	--output $output_file \
	--output-type u \
	--threads $threads \
	$genotypes_reconstituted_file \
	> $log_file

# Index the BCF file
bcftools index \
	--min-shift 30 \
	--threads $threads \
	$output_file

# Report to stdout the output file name and path.
echo -e "${GREEN}Annotation complete... Output file is located at : ${YELLOW}$output_file${NC}"
