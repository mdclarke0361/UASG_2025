#
rm(list=ls())
pdf(NULL)

#
suppressPackageStartupMessages({
	library(here)
	library(tidyverse)
    library(conflicted)
})

conflict_prefer_all("dplyr", quiet = TRUE)
conflicts_prefer(base::setdiff)

# Initialize the R env, get path variables for output.
source("01_source/initialize_script.r")

# Get logfile name
log_file <- assign_log_filename()

# Redirect output to the file and also keep it on the console
sink(file = log_file, split = TRUE)

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
clinical_data_file <- here(args[1])
genotype_report_file <- here(args[2])
expression_data_file <- here(args[3])

# Read in data
clinical_data <- read_csv(
    file = clinical_data_file
)

genotype_data <- read_tsv(
	file = genotype_report_file
)

expression_data <- read_rds(
	file = expression_data_file
)

# Construct a sample list for samples with complete data across clinical (mRS), RNA, and SNPs.
# Extract sample names from clinical data
clinical_samples <- clinical_data |>
    pull(sample_name)

# Find union between clinical data and SNP data.
genotype_samples <- genotype_data |>
    select(
        any_of(clinical_samples)
    ) |>
    colnames()
#

clinical_genotype_intersect <- intersect(clinical_samples, genotype_samples)

# Find intersection between all 3 datasets.
expression_samples <- expression_data |>
    select(
        !c(PROBEID, SYMBOL, GENENAME)
    ) |>
    colnames()

# Create a final sample list based on intersection of all 3 datasets
sample_list <- intersect(clinical_genotype_intersect, expression_samples)

#Filter datasets using intersection list.
clinical_data_filtered <- clinical_data |>
    filter(
        sample_name %in% sample_list
    )

genotype_data_filtered <- genotype_data |>
    select(
        1:8,
        all_of(sample_list)
    )

expression_data_filtered <- expression_data |>
    select(
        1:3,
        all_of(sample_list)
    )

# Pivot genotype data to combine it with clinical data and add genotype values under each rsID.
genotype_data_pivoted <- genotype_data_filtered |>
    select(
        all_of(sample_list),
        ID
    ) |>
    pivot_longer(
        all_of(sample_list),
        names_to = "sample_name",
        values_to = "alt_allele_count"
    ) |>
    pivot_wider(
        names_from = ID,
        values_from = alt_allele_count
    )

clinical_genotype_combined <- clinical_data_filtered |>
    left_join(
        genotype_data_pivoted
    )

# Save results as csv files
write_csv(
    clinical_genotype_combined,
    file = here(processed_data_dir, "patient_data.csv")
)

write_csv(
    expression_data_filtered,
    file = here(processed_data_dir, "expression_data.csv")
)

# Stop redirecting output
sink()