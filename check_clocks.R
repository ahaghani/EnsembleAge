#!/usr/bin/env Rscript

# Check what clocks we have
epiclocks <- readRDS("data/Clock coefficients.RDS")

cat("Clock families available:\n")
for (i in 1:length(epiclocks)) {
  cat(i, ".", names(epiclocks)[i], "\n")
  family_clocks <- epiclocks[[i]]
  for (j in 1:length(family_clocks)) {
    cat("  ", names(family_clocks)[j], "\n")
  }
  cat("\n")
}

# Check if we have the ensemble dynamic clocks
cat("Looking for ensemble dynamic clocks...\n")
found_dynamic <- FALSE
for (family_name in names(epiclocks)) {
  if (grepl("Ensemble.*Dynamic", family_name) || grepl("Dynamic", family_name)) {
    cat("Found potential dynamic clocks in:", family_name, "\n")
    found_dynamic <- TRUE
  }
}

if (!found_dynamic) {
  cat("No dynamic ensemble clocks found in coefficient file.\n")
  cat("Need to add EnsembleAge.Dynamic.Clocks.csv data to the package.\n")
}
