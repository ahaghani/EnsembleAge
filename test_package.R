#!/usr/bin/env Rscript

# Test script for EnsembleAge package
# This script tests the package functions with real data

# Load required libraries
library(dplyr)
library(tibble)
library(tidyr)
library(data.table)
library(plyr)
library(stringr)

# Source our package functions since we're testing before installation
source("R/transformation_functions.R")
source("R/clock3_functions.R")
source("R/platform_utils.R")
source("R/data_preprocessing.R")
source("R/predict_age.R")
source("R/examples.R")

cat("=== Testing EnsembleAge Package ===\n\n")

# Test 1: Mammal40k data
cat("Test 1: Loading Mammal40k data...\n")
mammal40k_path <- "/Users/ahaghani/Library/CloudStorage/GoogleDrive-ahaghani@altoslabs.com/My Drive/Amin documents/Steve projects/Research projects/Collaborative projects/Altos projects/Altos Data/N169.ET0340.AltosReddyCyno.m40k"

# Load normalized data
dat40k <- readRDS(file.path(mammal40k_path, "NormalizedData/all_probes_sesame_normalized.RDS"))
samps40k <- read.csv(file.path(mammal40k_path, "SampleSheetOutput_qy_final.csv"))

cat("Data dimensions:", nrow(dat40k), "probes x", ncol(dat40k) - 1, "samples\n")
cat("Sample sheet dimensions:", nrow(samps40k), "samples x", ncol(samps40k), "variables\n")

# Test platform detection
platform <- detect_platform(dat40k)
cat("Detected platform:", platform, "\n")

# Test data validation
validation <- validate_data_format(dat40k, samps40k)
cat("Data validation:", ifelse(validation$valid, "PASSED", "FAILED"), "\n")
if (!validation$valid) {
  cat("Issues found:\n")
  for (issue in validation$issues) {
    cat(" -", issue, "\n")
  }
}

# Test probe coverage
cat("Checking probe coverage...\n")
coverage <- check_probe_coverage(dat40k)
cat("Number of clocks with >50% coverage:", sum(coverage$Coverage_Percent > 50), "\n")
cat("Number of clocks with >80% coverage:", sum(coverage$Coverage_Percent > 80), "\n")

# Show top clocks by coverage
cat("Top 10 clocks by probe coverage:\n")
print(head(coverage, 10))

# Test preprocessing
cat("\nTesting preprocessing...\n")
processed <- preprocess_methylation_data(dat40k, samps40k, verbose = FALSE)
cat("Processed data:", nrow(processed$data), "probes x", nrow(processed$samples), "samples\n")

# Test age prediction
cat("\nTesting age prediction...\n")
tryCatch({
  results <- predictAgeAndAgeAcc(processed$data, processed$samples)
  cat("Prediction successful! Output dimensions:", nrow(results), "x", ncol(results), "\n")
  cat("Available predictions:\n")
  prediction_cols <- names(results)[grepl("epiAge", names(results))]
  for (i in 1:min(10, length(prediction_cols))) {
    cat(" -", prediction_cols[i], "\n")
  }
  
  # Save test results
  write.csv(results, "test_results_mammal40k.csv", row.names = FALSE)
  cat("Test results saved to test_results_mammal40k.csv\n")
  
}, error = function(e) {
  cat("Error in age prediction:", e$message, "\n")
})

cat("\n=== Test 1 Complete ===\n\n")

# Test 2: Simple wrapper function
cat("Test 2: Testing simple prediction wrapper...\n")
tryCatch({
  simple_results <- predict_age_simple(dat40k, samps40k, verbose = FALSE)
  cat("Simple prediction successful! Output dimensions:", nrow(simple_results), "x", ncol(simple_results), "\n")
}, error = function(e) {
  cat("Error in simple prediction:", e$message, "\n")
})

cat("\n=== All Tests Complete ===\n")
