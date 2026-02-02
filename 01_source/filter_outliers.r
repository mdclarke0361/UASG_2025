#
rm(list=ls())
pdf(NULL)

#
suppressPackageStartupMessages({
	library(here)
	library(tidyverse)
	library(conflicted)
	library(oligo)
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
raw_expression_data_file <- here(args[1])
outlier_list_file <- here(args[2])

# Read in the raw expression data from saved R object
raw_expression <- read_rds(
    raw_expression_data_file
)

# Get list of outliers from saved text file
outlier_list <- read_lines(
    outlier_list_file
)

# Filter eset object to exclude outliers
filtered_expression <- raw_expression[, !colnames(raw_expression) %in% outlier_list]

# In case of duplicated samples, deduplicate sample names, then keep only unique values.
deduplicated <- c()
samples <- sampleNames(filtered_expression)

for (sample in samples) {

  if (str_count(sample, "_") == 2) {

    sample <- str_remove(sample, "_[^_]*$")

  }

  deduplicated <- c(deduplicated, sample)

}

unique_samples <- unique(deduplicated)

# Subset the eset object to only include the deduplicated sample set.
deduplicated_expression <- filtered_expression[, sampleNames(filtered_expression) %in% unique_samples]

write_rds(
    deduplicated_expression,
    file = here(data_out, "filtered_expression_data.rds")
)

cat(
    "\nExpression data filtering complete. Filtered data saves as:",
    here(data_out, "filtered_expression_data.rds"),
    "\n"
)

# Stop redirecting output
sink()
