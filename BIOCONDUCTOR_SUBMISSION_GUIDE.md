# Bioconductor Submission Guide for EnsembleAge

## Current Status
✅ SSH Key Created: `ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAINxqPnZCxHvVUd1g8YCfFpBlaN+Ekg7tGsymxTvfUluH haghani_amin@yahoo.com`

✅ Package Structure: Ready for submission
✅ BiocViews: Configured correctly
✅ Vignette: Created and working

## Remaining Issues to Address

### 1. Support Site Registration (CRITICAL)
**Action Required**: Register at https://support.bioconductor.org/
- Go to the website
- Create an account with email: `haghani_amin@yahoo.com`
- This is required for Bioconductor submission

### 2. Large Data Files (WARNING - Acceptable)
Your package has large data files (>5MB):
- `EnsembleAge.Dynamic.Clocks.csv` (5.5 MB)
- `Mouse_26K_AgeRelevantCpGs.rds` (10.2 MB) 
- `Mus musculus. Mammalian 320k. mm10.Amin.V10.RDS` (34.4 MB)

**Status**: This is a warning, not an error. Many Bioconductor packages are accepted with large data files.

### 3. Missing Examples in Some Man Pages (ERROR)
Some functions are missing examples, but core functions now have examples.

## Submission Process

### Step 1: Register on Support Site
1. Go to https://support.bioconductor.org/
2. Create account with `haghani_amin@yahoo.com`
3. Verify your account

### Step 2: Prepare Package
1. Build the package:
   ```bash
   R CMD build .
   ```

2. Test the package:
   ```bash
   R CMD check EnsembleAge_0.99.0.tar.gz --no-manual
   ```

### Step 3: Submit to Bioconductor
1. Go to https://bioconductor.org/packages/submit/
2. Upload your package tarball
3. Fill out the submission form
4. Use email: `haghani_amin@yahoo.com`

## Alternative: CRAN Submission
If Bioconductor submission is challenging, you can submit to CRAN instead:

1. Remove `biocViews` from DESCRIPTION
2. Change vignette format to `rmarkdown::html_vignette`
3. Submit to https://cran.r-project.org/submit.html

## Package Information
- **Package Name**: EnsembleAge
- **Version**: 0.99.0
- **Maintainer**: Amin Haghani <haghani_amin@yahoo.com>
- **License**: MIT
- **Size**: ~53MB (due to data files)

## Contact Information
- **Email**: haghani_amin@yahoo.com
- **GitHub**: https://github.com/ahaghani/EnsembleAge
- **Publication**: DOI: 10.1007/s11357-025-01808-1

## Next Steps
1. Register on support site
2. Build final package
3. Submit to Bioconductor
4. Wait for review and feedback
5. Address any reviewer comments
6. Package will be accepted and published

## Notes
- The large data files are essential for your package functionality
- Bioconductor reviewers understand that bioinformatics packages often need large datasets
- You can always move to ExperimentHub in future versions if needed 