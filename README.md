# EnsembleAge

An R package for predicting epigenetic age using ensemble clock methods.

## Description

EnsembleAge provides tools for predicting epigenetic age using various clock methods including Universal Clocks, Elastic Epigenetic clocks, static ensemble methods, and dynamic ensemble clocks. The package is designed to work with methylation data and provides accurate age prediction across various tissues and species.

**This package accompanies the open access publication:**
> Haghani, A., et al. (2025). EnsembleAge: enhancing epigenetic age assessment with a multi‑clock framework. *GeroScience*. DOI: [10.1007/s11357-025-01808-1](https://doi.org/10.1007/s11357-025-01808-1)

The clocks and methods implemented in this package are freely available for research use under the open access license.

## Features

- Multiple epigenetic clock implementations
- Universal Clock 2 and 3 methods
- Elastic Epigenetic clocks with age transformation
- Static ensemble methods
- Support for multiple tissues and species
- Age acceleration calculations

## Installation

You can install the development version of EnsembleAge from GitHub:

```r
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Install EnsembleAge
devtools::install_github("ahaghani/EnsembleAge")
```

## Dependencies

The package requires the following R packages:
- dplyr
- tibble
- tidyr
- plyr
- stringr
- data.table

## Usage

### Quick Start

```r
library(EnsembleAge)

# Method 1: Simple prediction with automatic preprocessing
results <- predict_age_simple(your_methylation_data, your_sample_info)

# Method 2: Manual control over preprocessing
processed <- preprocess_methylation_data(your_methylation_data, your_sample_info)
results <- predictAgeAndAgeAcc(processed$data, processed$samples)
```

### Complete Example with Synthetic Data

```r
library(EnsembleAge)

# Create example data
example_data <- create_example_methylation_data(n_probes = 2000, n_samples = 30)
sample_info <- create_example_sample_info(sample_names = names(example_data)[-1])

# Run complete analysis
results <- run_example_analysis(n_probes = 2000, n_samples = 30)

# View results
head(results$age_predictions)
```

### Working with Different Platforms

```r
# Detect your platform
platform <- detect_platform(your_data)
cat("Detected platform:", platform)

# Check probe coverage
coverage <- check_probe_coverage(your_data)
print(coverage[coverage$Coverage_Percent > 80, ])  # Show well-covered clocks

# Validate data format
validation <- validate_data_format(your_data, your_samples)
if (!validation$valid) {
  print(validation$issues)
}
```

### Custom Preprocessing

```r
# Preprocess with custom settings
processed <- preprocess_methylation_data(
  dat0sesame = your_data,
  samps = your_samples,
  min_coverage = 0.7,           # Lower coverage threshold
  handle_missing = "impute",    # Impute missing values
  verbose = TRUE
)

# Run predictions
results <- predictAgeAndAgeAcc(processed$data, processed$samples)
```

## Data Requirements

### Methylation Data (`dat0sesame`)
- Must contain a `CGid` column with CpG site identifiers
- Sample columns containing methylation beta values (0-1)

### Sample Information (`samps`)
- `Basename`: Sample identifiers (required)
- `Age`: Chronological age (optional, defaults to 0)
- `SpeciesLatinName`: Species name (optional, defaults to "Mus musculus")
- `Female`: Sex indicator (optional)
- `Tissue`: Tissue type (optional)

## Clock Types

The package includes several types of epigenetic clocks:

1. **Universal Clocks**: Cross-species clocks (UniClock2, UniClock3)
2. **Elastic Epigenetic Clocks**: Age-transformed clocks
3. **Static Ensemble Clocks**: Ensemble methods combining multiple clocks
4. **Tissue-specific Clocks**: Specialized for different tissue types

## Output

The function returns a data frame with:
- Age predictions for various clock types
- Age acceleration values
- Clock-specific results for different tissues

## Author

Amin Haghani

## License

MIT

## Citation

If you use this package in your research, please cite:

```
Haghani, A., et al. (2025). EnsembleAge: enhancing epigenetic age assessment with a multi‑clock framework. GeroScience. 
DOI: 10.1007/s11357-025-01808-1
```

**BibTeX:**
```bibtex
@article{haghani2025ensembleage,
  title={EnsembleAge: enhancing epigenetic age assessment with a multi‑clock framework},
  author={Haghani, Amin and [other authors]},
  journal={GeroScience},
  year={2025},
  doi={10.1007/s11357-025-01808-1},
  url={https://doi.org/10.1007/s11357-025-01808-1}
}
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Issues

If you encounter any issues, please file them on the [GitHub issues page](https://github.com/ahaghani/EnsembleAge/issues).
