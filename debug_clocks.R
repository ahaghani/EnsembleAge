#!/usr/bin/env Rscript

# Debug script to see what's causing the duplicate keys issue
library(dplyr)
library(data.table)

# Load clock coefficients and check for potential duplicates after cleaning
epiclocks <- readRDS("data/Clock coefficients.RDS")

cat("Checking EnsembleAge.Dynamic clock names before and after cleaning:\n")
dynamic_clocks <- names(epiclocks[["EnsembleAge.Dynamic"]])

cat("\nOriginal names:\n")
for (i in 1:length(dynamic_clocks)) {
  cat(i, ":", dynamic_clocks[i], "\n")
}

cat("\nCleaned names:\n")
cleaned_names <- ifelse(grepl("^X[0-9]+\\.", dynamic_clocks), 
                       gsub("^X[0-9]+\\.", "", dynamic_clocks), 
                       dynamic_clocks)

for (i in 1:length(cleaned_names)) {
  cat(i, ":", cleaned_names[i], "\n")
}

cat("\nDuplicate cleaned names:\n")
duplicates <- cleaned_names[duplicated(cleaned_names)]
if (length(duplicates) > 0) {
  for (dup in unique(duplicates)) {
    cat("Duplicate:", dup, "\n")
    indices <- which(cleaned_names == dup)
    cat("  Original names:", paste(dynamic_clocks[indices], collapse = ", "), "\n")
  }
} else {
  cat("No duplicates found.\n")
}
