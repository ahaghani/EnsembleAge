# Script to help set up ExperimentHub for EnsembleAge package
# This script will help you prepare large data files for Bioconductor submission

library(ExperimentHub)
library(AnnotationHub)

# List of large files that need to be moved to ExperimentHub
large_files <- c(
  "data/EnsembleAge.Dynamic.Clocks.csv",
  "data/Mouse_26K_AgeRelevantCpGs.rds", 
  "data/Mus musculus. Mammalian 320k. mm10.Amin.V10.RDS"
)

# Check file sizes
cat("Large files that need to be moved to ExperimentHub:\n")
for(file in large_files) {
  if(file.exists(file)) {
    size_mb <- file.size(file) / (1024 * 1024)
    cat(sprintf("%s: %.1f MB\n", file, size_mb))
  }
}

cat("\nTo submit to Bioconductor, you have two options:\n\n")

cat("OPTION 1: Move data to ExperimentHub (Recommended)\n")
cat("1. Create an ExperimentHub package\n")
cat("2. Upload your data files to ExperimentHub\n")
cat("3. Update your package to load data from ExperimentHub\n")
cat("4. Remove large files from your package\n\n")

cat("OPTION 2: Submit to CRAN instead\n")
cat("1. Keep data files in package (CRAN allows larger packages)\n")
cat("2. Remove biocViews from DESCRIPTION\n")
cat("3. Submit to CRAN\n\n")

cat("For ExperimentHub setup, you'll need to:\n")
cat("1. Create a separate package for your data\n")
cat("2. Follow Bioconductor's ExperimentHub guidelines\n")
cat("3. Update your main package to load from ExperimentHub\n\n")

cat("Would you like me to help you with either option?\n") 