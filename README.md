# EnsembleAge

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D3.5.0-blue.svg)](https://www.r-project.org/)

**Multi-platform epigenetic age prediction using ensemble clock methods**

## 📖 Description

EnsembleAge provides tools for predicting epigenetic age using various clock methods including Universal Clocks, Elastic Epigenetic clocks, static ensemble methods, and dynamic ensemble clocks. The package **automatically detects your data platform** and prepares the data accordingly for accurate age prediction across various tissues and species.

**📄 This package accompanies the open access publication:**
> Haghani, A., et al. (2025). EnsembleAge: enhancing epigenetic age assessment with a multi‑clock framework. *GeroScience*. DOI: [10.1007/s11357-025-01808-1](https://doi.org/10.1007/s11357-025-01808-1)

The clocks and methods implemented in this package are freely available for research use under the MIT license.

## 🚀 Quick Start

```r
library(EnsembleAge)

# Basic usage - works with any supported platform automatically!
results <- predict_all_clocks(your_methylation_data, your_sample_sheet)

# View results by clock family
head(results)
```

## 📦 Installation

```r
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Install EnsembleAge
devtools::install_github("ahaghani/EnsembleAge")
```

### Dependencies

The package automatically installs these required R packages:
- `dplyr`, `tibble`, `tidyr`, `plyr`, `stringr`, `data.table`

## 🔧 Platform Support

EnsembleAge **automatically detects and handles** multiple methylation array platforms:

### 🧬 Supported Platforms

| Platform | Array Type | Probe Count | Auto-Detection | Status |
|----------|------------|-------------|----------------|--------|
| **Mammal320k** | Mammalian Methylation Array | ~320,000 | ✅ Automatic | ✅ Fully Supported |
| **Mammal40k** | Reduced Mammalian Array | ~40,000 | ✅ Automatic | ✅ Fully Supported |
| **Human EPIC** | Illumina EPIC Array | ~850,000 | ✅ Automatic | ✅ Fully Supported |
| **Human 450k** | Illumina 450k Array | ~450,000 | ✅ Automatic | ✅ Fully Supported |

### 🤖 What Happens Automatically

**For Human Data (EPIC/450k):**
1. ✅ Detects `Probe_ID` format and maps to `CGid` using EPIC annotation
2. ✅ Handles large datasets efficiently (filters ~850k → ~30k relevant probes)
3. ✅ Applies appropriate clocks: Universal, EnsembleDual (human-optimized)

**For Mammal320k Data:**
1. ✅ Detects probe suffix format (`cg123_BC21`) and maps to clean CpG IDs (`cg123`)
2. ✅ Uses Mammalian 320k annotation for accurate mapping
3. ✅ Handles duplicates by taking mean values

**For Mammal40k Data:**
1. ✅ Works directly with clean CpG IDs (`cg123456`)
2. ✅ Applies all mouse clocks: Static, Dynamic, Universal

**For All Platforms:**
- 🔄 **Orientation Detection**: Automatically detects if probes are rows vs columns
- 🧬 **Missing Probe Imputation**: Missing probes filled with 0.5 (neutral methylation)
- ⚡ **Efficient Loading**: Only loads probes required for clock predictions
- 📊 **Sample Matching**: Automatically matches samples between data and sample sheet

## 📋 Data Requirements

### Input Data Formats

**Methylation Data:**
- **Human**: Either `Probe_ID` (will auto-map) or `CGid` column
- **Mouse/Mammal**: Either `CGid` or platform-specific probe IDs (will auto-convert)
- **Values**: Beta values (0-1 range) for each sample

**Sample Sheet (CSV):**
- **Required**: `Basename` (sample names), `Age` (chronological age)
- **Optional**: `Female` (0/1), `SpeciesLatinName`, `Tissue`
- **Auto-filled**: Missing columns get sensible defaults with user notification

### 📝 Sample Sheet Template

Use our provided template:

```r
# Copy sample sheet template
sample_template <- read.csv(system.file("extdata", "sample_sheet_template.csv", package = "EnsembleAge"))
```

| Basename | Age | Female | SpeciesLatinName | Tissue |
|----------|-----|--------|------------------|--------|
| Sample1 | 25 | 0 | Mus musculus | Liver |
| Sample2 | 30 | 1 | Homo sapiens | Blood |

## 🎯 Clock Types & Usage

### 🔍 Specific Clock Predictions

```r
# Best for mouse data - main EnsembleAge clocks
mouse_results <- predict_ensemble_static(methylation_data, sample_sheet)

# Individual mouse clocks - 50 specialized clocks
dynamic_results <- predict_ensemble_dynamic(methylation_data, sample_sheet)

# Cross-species clocks (Universal + others)
universal_results <- predict_original_clocks(methylation_data, sample_sheet)

# Best for human data - dual-species optimized
human_results <- predict_ensemble_dual_static(methylation_data, sample_sheet)

# Everything - all available clocks
all_results <- predict_all_clocks(methylation_data, sample_sheet)
```

### 📊 Clock Families Included

| Clock Family | Best For | # Clocks | Description |
|--------------|----------|----------|-------------|
| **EnsembleAge.Static** | 🐭 Mouse (primary) | 2 | Main mouse ensemble clocks |
| **EnsembleAge.Dynamic** | 🐭 Mouse (detailed) | 50 | Individual specialized mouse clocks |
| **EnsembleDualAge.Static** | 👨 Human (primary) | 1 | Human-mouse optimized clock |
| **UniClock2/3** | 🌐 Cross-species | 6 | Universal mammalian clocks |
| **LifespanUberClock** | 🐭 Mouse variants | 12 | Lifespan-focused clocks |
| **DNAmAge*** | 🐭 Mouse categories | 39 | Development, Elastic, Intervention clocks |

## 💡 Usage Examples

### Example 1: Human Data (Automatic)

```r
library(EnsembleAge)

# Load your human methylation data (with Probe_ID or CGid)
# Package automatically detects format and maps appropriately
results <- predict_ensemble_dual_static(human_betas, human_samples)

# View human-optimized age predictions
head(results)
#   Basename Age epiClock                    epiAge AgeAccelation Female Tissue
# 1 Sample1   75 panTissue                  68.2      -6.8          0   Blood
# 2 Sample2   75 panTissue                  73.1      -1.9          1   Blood
```

### Example 2: Mouse Data (Automatic)

```r
# Load your mouse methylation data
# Package detects Mammal320k/40k format automatically
mouse_results <- predict_ensemble_static(mouse_betas, mouse_samples)

# Check what platform was detected
cat("Platform detected:", attr(mouse_results, "platform"))

# View main mouse age predictions
head(mouse_results[c("Basename", "Age", "epiAge", "AgeAccelation")])
```

### Example 3: All Clocks with Platform Info

```r
# Run all available clocks
all_results <- predict_all_clocks(methylation_data, sample_sheet, verbose = TRUE)

# See what clocks worked for your data
family_summary <- all_results %>%
  group_by(clockFamily) %>%
  summarise(
    n_clocks = length(unique(epiClock)),
    n_predictions = n(),
    mean_age = round(mean(epiAge, na.rm = TRUE), 1)
  )

print(family_summary)
```

### Example 4: Check Platform Compatibility

```r
# Check your data before running predictions
platform <- detect_platform(your_data)
cat("Detected platform:", platform)

# Check probe coverage for different clocks
coverage <- check_probe_coverage(your_data)
print(coverage[coverage$Coverage_Percent > 70, ])  # Well-covered clocks

# Validate your sample sheet
processed_samples <- prepare_sample_sheet(your_samples)
```

## 🔧 Advanced Options

### Custom Preprocessing

```r
# Manual control over preprocessing
processed <- preprocess_methylation_data(
  dat0sesame = your_data,
  samps = your_samples,
  min_coverage = 0.8,           # Coverage threshold
  handle_missing = "impute",    # Missing value strategy
  verbose = TRUE
)

# Run predictions on preprocessed data
results <- predictAgeAndAgeAcc(processed$data, processed$samples)
```

### Efficient Loading for Large Datasets

```r
# Recommended for large datasets (like EPIC arrays)
results <- predict_all_clocks(
  your_data, 
  your_samples,
  efficient_loading = TRUE,  # Only loads required probes
  verbose = TRUE
)
```

## 📈 Understanding Results

### Output Format

All prediction functions return a data frame with:

- **Basename**: Sample identifier
- **Age**: Input chronological age  
- **epiClock**: Clock name (e.g., "panTissue", "Blood", "Liver")
- **epiAge**: Predicted epigenetic age
- **AgeAccelation**: Age acceleration (epiAge - chronological age, residualized)
- **Female**: Sex (0=Male, 1=Female, NA=Unknown)
- **Tissue**: Tissue type
- **clockFamily**: Clock family name

### Choosing the Right Clocks

**For Human Samples:**
- 🥇 **Primary**: `predict_ensemble_dual_static()` - Best human accuracy
- 🥈 **Secondary**: `predict_original_clocks()` - Universal clocks

**For Mouse Samples:**
- 🥇 **Primary**: `predict_ensemble_static()` - Main mouse clocks  
- 🥈 **Detailed**: `predict_ensemble_dynamic()` - 50 specialized clocks
- 🥉 **Universal**: `predict_original_clocks()` - Cross-species clocks

**For All Data:**
- 📊 **Comprehensive**: `predict_all_clocks()` - Everything available

## 🛠️ Troubleshooting

### Common Issues

**"Cannot identify CpG ID column"**
- Solution: Ensure your data has either `CGid` or `Probe_ID` column

**"No matching samples found"**
- Solution: Check that sample names in data match `Basename` in sample sheet

**"Low probe coverage warnings"**
- Solution: This is normal; package automatically handles missing probes

**Platform detected incorrectly**
- Solution: Contact us with details; detection is usually accurate

## 👥 Author & Citation

**Author:** Amin Haghani (Altos Labs)  
**Email:** ahaghani@altoslabs.com

### Citation

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

## 📄 License

MIT License - Free for research and commercial use.

## 🤝 Contributing & Issues

- **Contributing**: Contributions are welcome! Please submit a Pull Request.
- **Issues**: Report bugs or request features on our [GitHub issues page](https://github.com/ahaghani/EnsembleAge/issues).
- **Questions**: Contact Amin Haghani at ahaghani@altoslabs.com

---

**🎯 TL;DR: Just run `predict_all_clocks(your_data, your_samples)` - everything else is automatic!**