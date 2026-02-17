# Initialize project scripts.

library(this.path)

# Ensure all project directories exist, if not, create them.

subdirectories <- c(
	"02_data/metadata", "02_data/processed", "02_data/raw", "02_data/reference",
	"03_results/figures", "03_results/reports", "03_results/tables", "03_results/logs"
)

for (dir in subdirectories) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# Get path for project directory
script_dir <- this.dir()
project_dir <- dirname(script_dir)

# Set variables to match repo file directory.
metadata_dir <- "02_data/metadata"
processed_data_dir <- "02_data/processed"
raw_data_dir <- "02_data/raw"
reference_data_dir <- "02_data/reference"
figure_out <- "03_results/figures"
report_out <- "03_results/reports"
table_out <- "03_results/tables"
log_out <- "03_results/logs"

# Assign log filename
assign_log_filename <- function() {
	script_name <- basename(this.path())
	script_basename <- sub("\\.[^.]*$", "", script_name)
	log_file <- file.path(log_out, paste0(script_basename, ".log"))
	return(log_file)
}

#
plot_colours <- c(
    "#4A6990", "#A73030", "#79AF97",
    "#DF8F44", "#6A6599", "#374E55",
    "#B1746F", "#8A8B79", "#7AA6DC",
    "#616530", "#642822", "#9A5324",
    "#0B1948"
)
