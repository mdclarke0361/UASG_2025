# UASG 2025

# 1. Set up

## Clone Repo

Clone this repo to desired project directory.

## Recreate Conda Env

Below is an example command. Create conda environment from environment.yml in main project directory.
Use desired project name in place of env_name.
```bash
env_name="UASG"

conda env create \
	--name $env_name \
	--file environment.yml
```

All scripts reference paths based on the directory of that script. It is important that all scripts remain in the 01_source directory.

# 2. Processing SNP Data

If required, use Genome Studio to create a custom cluster file: https://www.illumina.com/content/dam/illumina/gcs/assembled-assets/marketing-literature/custom-cluster-file-tech-note-m-gl-02142/custom-cluster-file-tech-note-m-gl-02142.pdf

This project uses a custom array and has therefore required this step.

Conversion to GTC and VCF can be done using Illumina's DRAGEN or a GATK tool (Broad Institute). DRAGEN was chosen for this workflow as it is currently the method recommended by Illumina. Download the DRAGEN software from Illumina website and Unzip the downloaded file and save into the /source directory.

DRAGEN also converts the GTP files to VCF for use with downstream processes. Note: this requires the .csv version of the manifest file (supplied by illumina along with the bpm manifest).

Requirements:
- Place all input data from Illumina, including IDAT files, manifest files and cluster file into the project's data directory.

The following files are project-specific, and have been placed in the project's data directory. They will be required by commands within this workflow.
- manifest file (.bpm)
- sample metadata (.csv)
- sample manifest (.csv)
- A directory of IDAT files ie. 02_data/raw_data/idat_files

For this part of the workflow, the steps are outlined as follows:
1. genotyping - get allele calls from the raw intensity data.
2. imputation - run imputation to get allele calls from missing markers.
3. annotation - add gene identifiers to each marker.
4. reporting - preparation of data for statistical analysis.

## 2a. Genotype Calling

Download the recommended fasta and index files from Illumina website. 
The below example is done for GRCh37 and placed within the project reference data directory:
```bash
data_ref_dir="02_data/reference"
dl_file="${data_ref_dir}/ref_genome.zip"

target_url="https://webdata.illumina.com/downloads/productfiles/microarray-analytics-array/GRCh37_genome.zip"

wget \
	-O $dl_file \
	$target_url

unzip -d  $data_ref_dir $dl_file &&
	rm $dl_file
```

Run the script from source directory to get genotype calls from raw data
```bash
# Assign file paths to dependecies
bpm_manifest="02_data/metadata/snp_array_manifest.bpm"
cluster_file="02_data/metadata/cluster_file.egt"
idat_dir="02_data/raw/idat_files"
sample_manifest="02_data/metadata/snp_sample_manifest.csv"
csv_manifest="02_data/metadata/snp_array_manifest.csv"
ref_genome_fasta="02_data/reference/GRCh37_genome.fa"

# Run genotyping script
01_source/genotyping.sh \
	$bpm_manifest \
	$cluster_file \
	$idat_dir \
	$sample_manifest \
	$csv_manifest \
	$ref_genome_fasta
```

In this project the file header assigns names to each sample based on the Illumina barcode for that sample run. The below code will reference the sample manifest file and replace these barcoded headers with sample names used in our study.

```bash
# Assign file paths to dependecies
genotypes_file="03_results/processed/genotypes.bcf"
sample_manifest="02_data/metadata/snp_sample_manifest.csv"

# Run reheader script
01_source/reheader.sh \
	$genotypes_file \
	$sample_manifest
```

Filtering out internal controls based on sample name, then normalizing the outputs with bcftools
In this study, control samples have the prefix 'NA'.

```bash
# Assign file paths to dependecies
genotypes_fixed_header_file="02_data/processed/genotypes_fixed_header.bcf"
sample_list_file="02_data/metadata/snp_sample_names.txt"
control_sample_prefix="NA"

# Run filtering script
01_source/filter_controls.sh \
	$genotypes_fixed_header_file \
	$sample_list_file \
	$control_sample_prefix
```

```bash
# Assign file paths to dependecies
genotypes_filtered_controls_file="02_data/processed/genotypes_filtered_controls.bcf"
ref_genome_fasta="02_data/reference/GRCh37_genome.fa"

# Run normalization script
01_source/normalize.sh \
	$genotypes_filtered_controls_file \
	$ref_genome_fasta
```

## 2b. Imputation

Phasing and Imputation of Missing Genotypes with BEAGLE 5.5

Notes:
- Only autosomal and X chromosomes are kept for the imputation.
- The %INFO field will be populated with the imputation results.

http://faculty.washington.edu/browning/beagle/beagle.html

Get reference VCF files and the accompanying index files from the BEAGLE website.
```bash
# Set directory for the reference vcf files to be saved
ref_vcf_dir="02_data/reference/ref_vcf_files"
target_url="https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/"

wget \
	-r \
	--no-parent \
	--reject "index.html*" \
	--no-check-certificate \
	--directory-prefix $ref_vcf_dir \
	--no-directories \
	$target_url
```

The reference vcf files are split by chromosome. To best incorporate them with this workflow, we can concatenate the vcf files and re-index the concatenated file.

We can test the chromosome names of our bead chip data against the chromosome names in the reference. They must match for the imputation software to run.

```bash
# Assign file paths to dependecies
genotypes_normalized_file="02_data/processed/genotypes_normalized.bcf"
ref_vcf_dir="02_data/reference/ref_vcf_files"

# Run script to compile reference vcf file
01_source/compile_reference.sh \
	$genotypes_normalized_file \
	$ref_vcf_dir
```

The lists of chromosome names from each file are logged in 03_results/processed/log_files.
If necessary, convert the chromosome names in either file. This can be done with BCFTools, but is not carried out in this workflow.

Separate the chromosomes that don't have a reference from the sample BCF file. Convert back to a compressed VCF for use with BEAGLE.

The genotype call file needs to be prepared for imputation by filtering out incompatable chromosomes and converting X chromosome markers to haploid. See the imputation software for more information on why this is necassary.

```bash
#Assign file paths to dependecies
genotypes_normalized_file="02_data/processed/genotypes_normalized.bcf"
chr_list_file="03_results/reports/genotypes_chr_list.txt"
chr_to_remove="Y,MT"

01_source/filter_chr.sh \
	$genotypes_normalized_file \
	$chr_list_file \
	$chr_to_remove
```

Convert Haploid X Chromosome Markers to Diploid
Correct for male haptoid genotypes on chromosome X. BEAGLE only accepts diploid values.

```bash
genotypes_filtered_chr_file="02_data/processed/genotypes_filtered_chr.vcf.gz"
chr_length="155260478"

01_source/correct_ploidy.sh \
	$genotypes_filtered_chr_file \
	$chr_length
```

Getting reference .map files required by BEAGLE 5.5 for imputation.
Files for hg19 (GRCh37) can be downloaded from the BEAGLE site.

Since the unzipped .map files do not contain headers, and are in text format, we can concatenated them with `cat`.
Note that we must sort correctly to maintain chromosome order.

Pseudoautosomal regions (Par_1 and Par_2) would need to be run separately in BEAGLE. They will be left out of this analysis.

```bash
#
map_file_dir="02_data/reference/plink_map_files"
target_url="https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip"

wget \
	-P $map_file_dir \
	--no-check-certificate \
	$target_url

unzip "${map_file_dir}/plink.GRCh37.map.zip" &&
rm "${map_file_dir}/plink.GRCh37.map.zip"

# Create new file for output
concatenated_map_file="02_data/reference/concatenated_map_file.map"
touch $concatenated_map_file

for i in {1..22} X; do

	cat "${map_file_dir}/plink.chr${i}.GRCh37.map" >> $concatenated_map_file

done
```

Use BEAGLE 5.5 to run phasing and imputation.
Beagle is installed through bioconda, executed from bash, but runs in java. Note that there is an additional argument - '-Xmx<Gb>b' to specify the total memory to allocate for this process.

Set a threshold for the quality of imputations to be retained. Filter the files based on DR2 value to retain only quality imputations. Convert back to BCF format.

We used a DR2 threshold of 0.8 for this study.

```bash
# Based on available memory, specify the number in GB to allocate.
total_available_mem=$(free -g | grep "Mem:" | awk '{print $7}')
memory_to_use=$(echo "${total_available_mem}*0.8" | bc -l) # Use 80% of available mem
memory_to_use=${memory_to_use%.*} # Round to whole number

#
genotypes_fixed_ploidy_file="02_data/processed/genotypes_fixed_ploidy.vcf.gz"
concatenated_ref_vcf="02_data/reference/concatenated_ref.vcf.gz"
concatenated_map_file="02_data/reference/concatenated_map_file.map"
dr2_threshold=0.8

#
01_source/impute.sh \
	$memory_to_use \
	$genotypes_fixed_ploidy_file \
	$concatenated_ref_vcf \
	$concatenated_map_file
```

The vcf file 'INFO' field will now contain:
'DR2' subfield with the estimated squared correlation between the estimated allele dose and the true allele dose.
'AF' subfield with the estimated alternate allele frequencies in the target samples.
'IMP' flag if the marker is imputed.

The Genotype is accompanied by the allele dose and genotype probabilities.
'DS' The sum of the two allele probabilities (the allele dose).
'GP' The probilities for each of 0|0, 1|0, and 1|1 genotypes.

BEAGLE only outputs markers that are present in both the reference VCF and the sample VCF. If measured markers have been omitted in the imputed dataset, then we can extract the imputed markers from the imputed dataset and combine them with the measured dataset to reconstitute the full dataset; including markers from the Y and MT chromosomes.

The last genotype data file prior to removing the Y and MT chromosomes was produced from the normalization step.

```bash
# 
genotypes_imputed_filtered_file="03_results/processed/genotypes_imputed_filtered.bcf"
genotypes_normalized_file="03_results/processed/genotypes_normalized.bcf"

#
01_source/reconstitute.sh \
	$genotypes_imputed_filtered_file \
	$genotypes_normalized_file
```

Non-imputed markers will have a FORMAT of GT:GS:BAF:LRR, while the imputed markers are GT:DS.

## 2c. Gene Annotation

BCFTools will allow us to annotate the data with gene names for each SNP.
Adding gene annotations as per discussion on: https://www.biostars.org/p/122690/

First, a gene list must be acquired to add to the data.

For use with BCFTools we will need to compress and index the gene lists. Indexing of BED files can be done with Tabix.
```bash
gene_list_file="02_data/reference/gene_list_GRCh37.bed"
sorted_gene_list_file="02_data/reference/sorted_gene_list_GRCh37.bed"
target_url="https://www.cog-genomics.org/static/bin/plink/glist-hg19"

wget \
	-O $gene_list_file \
	$target_url

# Gene list can be sorted by chromosome, then numerically sorted by start position.
# Change the space-delimited file to tab-delimited as per BED specification.
sort -k1,1 -k2,2n $gene_list_file | \
	tr ' ' '\t' >  $sorted_gene_list_file

# Clean up temp file.
rm $gene_list_file

# Create a compressed copy of the gene list for use with bcftools, then index the gene list.
bgzip \
	--keep \
	$sorted_gene_list_file

tabix "${sorted_gene_list_file}.gz"
```

Use BCFTools to annotate the genotype file.
Gene name will be written as 'GENE=xxx' in the INFO column.
```bash
#
genotypes_reconstituted_file="03_results/processed/genotypes_reconstituted.bcf"
gene_list_file="02_data/reference/gene_list_GRCh37.bed.gz"

#
01_source/annotate.sh \
	$genotypes_reconstituted_file \
	$gene_list_file
```

Rename the final processed file for easy recognition during analysis.
```bash
#
genotypes_annotated_file="03_results/processed/genotypes_annotated.bcf"
genotypes_annotated_index_file="03_results/processed/genotypes_annotated.bcf.csi"
genotypes_processed_file="03_results/processed/genotypes_processed.bcf"
genotypes_processed_index_file="03_results/processed/genotypes_processed.bcf.csi"

mv -n $genotypes_annotated_file $genotypes_processed_file
mv -n $genotypes_annotated_index_file $genotypes_processed_index_file
```

## 2d. Compile Reference SNP List

Download reference file from UCSC database for list of variants and positions. Convert the big bed format to a standard BED file.

The file comes from the online database as a 'BigBed' file and must be converted.

Correct the chromosome names in the reference file and extract only positions and rsIDs.
```bash
#
snp_list_bb_file="02_data/reference/dbsnp153common_temp.bb"
target_url="http://hgdownload.gi.ucsc.edu/gbdb/hg19/snp/dbSnp153Common.bb"
snp_list_bed_file="02_data/reference/dbsnp153common_temp.bed"
output_file="02_data/reference/dbsnp_153common.bed"

#
wget \
	-O $snp_list_bb_file
	$target_url

bigBedToBed $snp_list_bb_file $snp_list_bed_file

# Set flag to indicate file header
header=TRUE

# Loop through the bed file and edit line by line
while read -r line; do

	read -ra snp_line <<< "$(echo $line)"

	# Assign variables from array elements.
	chr=${snp_line[0]}
	start=${snp_line[1]}
	end=${snp_line[2]}
	rsID=${snp_line[3]}

	if [[ header -eq "TRUE" ]]; then

		echo $line | awk -v OFS='\t' '{print $1, $2, $3, $4}' > $output_file

		header="FALSE"

	else

		new_name=${chr//"chr"/""}

		echo -e "$new_name\t$start\t$end\t$rsID" >> $output_file
		
	fi

done < $snp_list_bed_file

# Clean up temporary files
rm $snp_list_bb_file $snp_list_bed_file
```

# 3. Processing RNA Data

Processing Affymetrix Clariom S Human Data for UASG

For this part of the workflow, the data processing is done in R. The steps are outlined as follows:
1. Read Raw Data - read the .cel files and perform initial quality control.
2. Filtering Outliers and Duplicates
3. RMA Normalization
4. Annotate Markers with Gene Names

## 3a. Process Raw Data

```bash
cel_file_dir="02_data/raw/cel_files"
rna_sample_manifest="02_data/metadata/rna_sample_manifest.csv"

Rscript 01_source/process_raw_data.r \
	$cel_file_dir \
	$rna_sample_manifest
```

## 3b. Filter Raw Data
A list of outliers is saved to the output directory and can then be used for filtering out the named samples as required.

We filter out all outliers prior to normalizing the data.

Note: This data set contains some duplicated samples. The following code calls an R script which will deduplicate the dataset. This may not be compatible with other datasets and should be removed from the R script in that case.

```bash
raw_expression_data_file="02_data/processed/raw_expression_data.rds"
outlier_list_file="03_results/reports/expression_outlier_list.txt"

Rscript 01_source/filter_outliers.r \
	$raw_expression_data_file \
	$outlier_list_file
```

## 3c. RMA Normalization

Robust Multichip Average preprocessing methodology. 
This strategy allows background subtraction, quantile normalization and summarization (via median-polish).

```bash
#
filtered_expression_data_file="02_data/processed/filtered_expression_data.rds"

Rscript 01_source/rma_normalization.r \
	$filtered_expression_data_file
```

## 3d. Annotation

Note: This will remove probes that do not have an assigned gene symbol.

```bash
#
normalized_expression_data_file="02_data/processed/normalized_expression_data.rds"

#
Rscript 01_source/annotate.r \
	$normalized_expression_data_file \
	$output_dir
```

# 4. Analysis 1: CD40 5077(T) Polymorphism (rs1883832)

Export genotype data on the target SNP
```bash
target_snp_list="02_data/reference/target_snp_list.txt"
processed_genotype_file="03_results/processed/genotypes_processed.bcf"
ref_snp_list="02_data/reference/dbsnp_153common.bed"
report_file="03_results/reports/analysis_1_snp_report.tsv"

01_source/analysis/target_snp_report.sh \
	$target_snp_list \
	$processed_genotype_file \
	$ref_snp_list \
	$report_file
```

Combine datasets and get a single sample set that has complete genotype, clinical, and expression data.

This assumes that the clinical data is a complete list of the samples to be studied. Therefore the sample list used as a reference is taken from the clinical data.
```bash
clinical_data_file="02_data/metadata/clinical_data.csv"
genotype_report_file="04_analysis/analysis_1/target_snp_report.tsv"
expression_data_file="03_results/processing_rna/annotated_expression_data.rds"

Rscript 01_source/analysis/combine_datasets.r \
	$clinical_data_file \
	$genotype_report_file \
	$expression_data_file \
```

Analysis is being completed in an interactive terminal and documented in an .rmd file:

## 4b. Export SNP Data on Target SNPs

Create a text file with rsIDs of the target SNPs. Add the rsID one per line.
```bash
#
analysis_name="analysis_01"
analysis_dir=04_analysis/${analysis_name}

mkdir -p $analysis_dir 2> /dev/null

#
target_snp_list=${analysis_dir}/target_snp_list.txt
touch $target_snp_list

#
processed_genotype_file=03_results/processed/genotypes_processed.bcf
ref_snp_list=02_data/reference/dbsnp_153common.bed

# Create target SNP report.
01_source/analysis/target_snp_report.sh \
	$analysis_dir \
	$target_snp_list \
	$processed_genotype_file \
	$ref_snp_list
```

## 4c. Export Data on Target Genes

Create a text file with gene symbols of the target genes. Add the gene symbol one per line.
```bash
#
analysis_name="analysis_01"
analysis_dir=04_analysis/${analysis_name}

mkdir -p $analysis_dir 2> /dev/null

#
target_gene_list=${analysis_dir}/target_gene_list.txt
touch $target_gene_list

#
processed_genotype_file=03_results/processed/genotypes_processed.bcf
ref_gene_list=02_data/reference/sorted_gene_list_GRCh37.bed

# Set buffer values for before and after locus for including markers near the gene.
bp_before_gene=250
bp_after_gene=250

# Create target gene report.
01_source/analysis/target_gene_report.sh \
	$analysis_dir \
	$target_gene_list \
	$processed_genotype_file \
	$ref_gene_list \
	$bp_before_gene \
	$bp_after_gene
```
