#!/usr/bin/env Rscript

# Script to integrate the ensemble clocks from CSV files into the main coefficient structure

library(dplyr)

# Load the existing clock coefficients
epiclocks <- readRDS("data/Clock coefficients.RDS")

# Function to convert CSV clock data to the coefficient format
csv_to_clock_format <- function(csv_file) {
  clock_data <- read.csv(csv_file, stringsAsFactors = FALSE)
  
  # Get all clock names (columns except CGid)
  clock_names <- names(clock_data)[names(clock_data) != "CGid"]
  
  clocks <- list()
  
  for (clock_name in clock_names) {
    # Get non-zero coefficients
    coef_data <- clock_data[, c("CGid", clock_name)]
    names(coef_data) <- c("CGid", "Coef")
    
    # Remove rows with zero coefficients (except Intercept)
    coef_data <- coef_data[coef_data$Coef != 0 | coef_data$CGid == "Intercept", ]
    
    clocks[[clock_name]] <- list(
      CGid = coef_data$CGid,
      Coef = coef_data$Coef
    )
  }
  
  return(clocks)
}

# Add EnsembleAge.Dynamic.Clocks
cat("Adding EnsembleAge.Dynamic.Clocks...\n")
dynamic_clocks <- csv_to_clock_format("data/EnsembleAge.Dynamic.Clocks.csv")
epiclocks[["EnsembleAge.Dynamic"]] <- dynamic_clocks

# Add EnsembleDualAge.Static.Clocks
cat("Adding EnsembleDualAge.Static.Clocks...\n")
dual_clocks <- csv_to_clock_format("data/EnsembleDualAge.Static.Clocks.csv")
epiclocks[["EnsembleDualAge.Static"]] <- dual_clocks

# Save the updated clock coefficients
saveRDS(epiclocks, "data/Clock coefficients.RDS")

cat("Updated clock coefficients saved!\n")
cat("New clock families added:\n")
cat("- EnsembleAge.Dynamic with", length(dynamic_clocks), "clocks\n")
cat("- EnsembleDualAge.Static with", length(dual_clocks), "clocks\n")

# Display the dynamic clock names
cat("\nEnsembleAge.Dynamic clocks:\n")
for (clock_name in names(dynamic_clocks)) {
  cat(" -", clock_name, "\n")
}

cat("\nEnsembleDualAge.Static clocks:\n")
for (clock_name in names(dual_clocks)) {
  cat(" -", clock_name, "\n")
}
