#
rm(list=ls())
pdf(NULL)

#
suppressPackageStartupMessages({
	library(here)
	library(tidyverse)
	library(conflicted)
	library(oligo)
    library(pd.clariom.s.human.ht)
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
filtered_expression_data_file <- here(args[1])

# Read in the raw expression data from saved R object
filtered_expression <- read_rds(
    filtered_expression_data_file
)

# Use RMA normalization
normalized_expression <- rma(filtered_expression)

# Get visuals of data distribution of normalized data.
hist <- oligo::hist(normalized_expression)

# Determine consistency of intensity values using oligo::boxplot
png(
	filename = here(figure_out, "normalized_intensity_boxplot.png"),
	width = 1000, 
	height = 1000
)

oligo::boxplot(
    normalized_expression, 
    main = "Normalized Expression Intensities", 
    targets = "core"
)

dev.off()

#
png(
	filename = here(figure_out, "normalized_intensity_density.png"),
	width = 1000, 
	height = 1000
)

oligo::hist(
	normalized_expression,
	main = "Normalized Expression Density",
)

dev.off()

# Run PCA to check dispersion and identify outliers

# Extract expression values and log-transform (essential for PCA)
expression_values <- exprs(normalized_expression)
log_data <- log2(expression_values + 1)

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

# Define outliers as any point > 2 Standard Deviations from the mean on PC1 or PC2
mean_pc1 <- mean(pca_df$PC1)
sd_pc1 <- sd(pca_df$PC1)
mean_pc2 <- mean(pca_df$PC2)
sd_pc2 <- sd(pca_df$PC2)

pca_df$is_outlier <- abs(pca_df$PC1 - mean_pc1) > 2*sd_pc1 | 
	abs(pca_df$PC2 - mean_pc2) > 2*sd_pc2

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
		xintercept = c(mean_pc1 + 2*sd_pc1, mean_pc1 - 2*sd_pc1), 
		linetype="dashed", 
		color="gray"
	) +
	geom_hline(
		yintercept = c(mean_pc2 + 2*sd_pc2, mean_pc2 - 2*sd_pc2), 
		linetype="dashed", 
		color="gray"
	) +
	# Aesthetics
	theme_minimal()

#
ggsave(
	plot,
	height = 10,
	width = 10,
	filename = here(figure_out, "normalized_expression_pca.png")
)

# Save normalized data as R object
write_rds(
	normalized_expression,
	file = here(processed_data_dir, "normalized_expression_data.rds")
)

cat(
	"Normalization complete. Data is saved as an R object at:",
	here(output_dir, "raw_expression_data.rds"),
	"\n"
)

# Stop redirecting output
sink()