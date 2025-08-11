# EnsembleAge

An R package for predicting epigenetic age using ensemble clock methods.

## Description

EnsembleAge provides tools for predicting epigenetic age using various clock methods including Universal Clocks, Elastic Epigenetic clocks, and static ensemble methods. The package is designed to work with methylation data and provides accurate age prediction across various tissues and species.

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
devtools::install_github("your-username/EnsembleAge")
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

```r
library(EnsembleAge)

# Load your methylation data
# dat0sesame should be a data frame with CGid column and sample columns
# samps should be a data frame with sample information

# Predict ages
results <- predictAgeAndAgeAcc(dat0sesame, samps)
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

[Citation information to be added]

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Issues

If you encounter any issues, please file them on the [GitHub issues page](https://github.com/your-username/EnsembleAge/issues).
