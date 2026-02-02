#
rm(list=ls())
pdf(NULL)

#
suppressPackageStartupMessages({
	library(here)
	library(tidyverse)
	library(conflicted)
	library(oligo)
  library(clariomshumantranscriptcluster.db)
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
normalized_expression_data_file <- here(args[1])

# Read in the expression data from saved R object
normalized_expression <- read_rds(
  normalized_expression_data_file
)

# Convert eset object to tibble
expression_data <- exprs(normalized_expression) |>
  as_tibble() |>
  # Create a column for probe IDs
  mutate(
    PROBEID = rownames(normalized_expression),
    .before = 1
  )

# Bring in database for gene symbols and names
anno_db <- clariomshumantranscriptcluster.db
probe_id_list <- featureNames(normalized_expression)

# Lookup gene symbols and names from reference
mapped_gene_list <- AnnotationDbi::select(
  anno_db, 
  keys = probe_id_list, 
  columns = c("SYMBOL", "GENENAME"), 
  keytype = "PROBEID"
  )

# Remove probes without gene symbols
filtered_gene_list <- mapped_gene_list |>
	filter(
		!is.na(mapped_gene_list$SYMBOL)
	) |>
  as_tibble()

#
annotated_expression <- filtered_gene_list |>
  left_join(
    expression_data
  )

# Save Data as R Objects for Further Processing Steps
write_rds(
  annotated_expression,
  file = here(processed_data_dir, "annotated_expression_data.rds")
)

cat(
	"Annotation complete. Annotated data is saved as an R object at:",
	here(processed_data_dir, "annotated_expression_data.rds"),
  "\n"
)

# Stop redirecting output
sink()
