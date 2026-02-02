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
cel_file_dir <- here(args[1])
rna_sample_manifest <- here(args[2])

# Load in CEL files and arrange to incorporate appropriate sample names.
# Assemble list of CEL files to be processed.
cel_file_list <- list.celfiles(cel_file_dir, full.names = TRUE)

# Load in sample manifest
sample_manifest <- read_csv(
  file = rna_sample_manifest
)

# Lookup CEL filename in manifest to create list of sample IDs.
manifest <- tibble(
  filename = character(),
  array = character()
  )

for (filename in cel_file_list) {
  array <- sub(".CEL", "", basename(filename))
  array <- sub("A_", "B_", array)
  manifest <- manifest |>
    add_row(
      filename = filename,
      array = array
    )
}

# Use tibble to only load in patient samples (omit controls).
# modify duplicated samples to increment second value
sample_file_df <- sample_manifest |>
  left_join(
    manifest
  ) |>
  filter(
    ! is.na(array)
  ) |>
  mutate(
    sample_id = make.unique(sample_id, sep = "_")
  )

# Load in raw data (CEL files)
raw_expression <- read.celfiles(
  sample_file_df$filename,
  sampleNames = sample_file_df$sample_id
)

#
hist <- oligo::hist(raw_expression)

# Determine consistency of intensity values using oligo::boxplot
png(
	filename = here(figure_out, "raw_intensity_boxplot.png"),
	width = 1000, 
	height = 1000
)

oligo::boxplot(
    raw_expression, 
    main = "Raw Expression Intensities", 
    targets = "core"
)

#
png(
	filename = here(figure_out, "raw_intensity_density.png"),
	width = 1000, 
	height = 1000
)

oligo::hist(
	raw_expression,
	main = "Raw Expression Density"
)

dev.off()

# Run PCA to check dispersion and identify outliers

# Extract raw expression values and log-transform (essential for PCA)
exp_raw <- exprs(raw_expression)
log_data <- log2(exp_raw + 1)

# Must transpose data to make samples rows.
transposed_log_data <- t(log_data)

# Calculate PCA (transpose data so samples are rows)
pca <- prcomp(
	transposed_log_data, 
	center = TRUE, 
	scale. = TRUE
)

# Convert PCA data into a dataframe for plotting.
pca_df <- as.data.frame(pca$x)
pca_df$sample <- rownames(pca_df)

# Set outlier threshold as x * standard deviations
oulier_threshold <- 2.5

# Define outliers
mean_pc1 <- mean(pca_df$PC1)
sd_pc1 <- sd(pca_df$PC1)
mean_pc2 <- mean(pca_df$PC2)
sd_pc2 <- sd(pca_df$PC2)

pca_df$is_outlier <- abs(pca_df$PC1 - mean_pc1) > oulier_threshold*sd_pc1 | 
	abs(pca_df$PC2 - mean_pc2) > oulier_threshold*sd_pc2


# Plot with ggplot
plot <- pca_df |>
	ggplot(
		aes(
			x=PC1,
			y=PC2,
			label=sample
		)
	) +
	# Draw the points, coloring them if they are outliers
	geom_point(
		aes(
			color=is_outlier
		), 
		size=3
	) +
	# Add labels (only to outliers to keep it clean, or remove 'data=...' to label all)
	geom_text(
		data=subset(pca_df, is_outlier==TRUE), 
		vjust=-0.5, 
		size=2
	) +
	# Add dashed lines to show the 2 SD threshold
	geom_vline(
		xintercept = c(mean_pc1 + oulier_threshold*sd_pc1, mean_pc1 - oulier_threshold*sd_pc1), 
		linetype="dashed", 
		color="gray"
	) +
	geom_hline(
		yintercept = c(mean_pc2 + oulier_threshold*sd_pc2, mean_pc2 - oulier_threshold*sd_pc2), 
		linetype="dashed", 
		color="gray"
	) +
	# Aesthetics
	theme_minimal()

plot

#
ggsave(
	plot,
	height = 10,
	width = 10,
	filename = here(figure_out, "raw_expression_pca.png")
)

# Save a list of outliers for reference and for filtering of data
outlier_list <- pca_df |>
	filter(
        is_outlier == TRUE
    ) |>
    pull(
        sample
    )

write(
	outlier_list,
	file = here(report_out, "expression_outlier_list.txt")
)

# Save raw data as R object
write_rds(
	raw_expression,
	file = here(data_out, "raw_expression_data.rds")
)

cat(
	"Processing complete. Outliers are listed in:",
	here(report_out, "outlier_list.txt"),
	"\nRaw data is saved as an R object at:",
	here(data_out, "raw_expression_data.rds"),
	"\n"
)

# Stop redirecting output
sink()