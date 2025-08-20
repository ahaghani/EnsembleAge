#' Predict Age from SummarizedExperiment Object
#'
#' This function demonstrates interoperability with Bioconductor's SummarizedExperiment
#' class, which is commonly used for storing methylation data from packages like
#' sesame and minfi.
#'
#' @param se A SummarizedExperiment object containing methylation beta values
#' @param sample_sheet Optional data frame with sample information. If NULL,
#'   will attempt to extract from colData(se)
#' @param assay_name Name of the assay containing beta values (default: "beta")
#' @param verbose Logical indicating whether to print progress messages
#' @return Data frame with age predictions
#' @export
#' @examples
#' \donttest{
#' # Create example SummarizedExperiment
#' library(SummarizedExperiment)
#' 
#' # Simulate methylation data
#' n_probes <- 100
#' n_samples <- 6
#' beta_values <- matrix(rbeta(n_probes * n_samples, 2, 2), 
#'                       nrow = n_probes, ncol = n_samples)
#' rownames(beta_values) <- paste0("cg", sprintf("%08d", 1:n_probes))
#' colnames(beta_values) <- paste0("Sample_", 1:n_samples)
#' 
#' # Create sample metadata
#' col_data <- data.frame(
#'   Basename = paste0("Sample_", 1:n_samples),
#'   Age = c(0.5, 1.2, 2.1, 0.8, 1.5, 2.8),
#'   Female = c(1, 0, 1, 1, 0, 1),
#'   Tissue = "Blood"
#' )
#' 
#' # Create SummarizedExperiment
#' se <- SummarizedExperiment(
#'   assays = list(beta = beta_values),
#'   colData = col_data
#' )
#' 
#' # Predict ages
#' results <- predict_age_from_se(se)
#' }
predict_age_from_se <- function(se, sample_sheet = NULL, assay_name = "beta", verbose = TRUE) {
  
  # Check if SummarizedExperiment is loaded
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is required but not installed.")
  }
  
  # Validate input
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment object")
  }
  
  if (verbose) {
    cat("Processing SummarizedExperiment object...\n")
    cat("Dimensions:", nrow(se), "features x", ncol(se), "samples\n")
  }
  
  # Extract beta values
  if (!assay_name %in% names(SummarizedExperiment::assays(se))) {
    stop("Assay '", assay_name, "' not found. Available assays: ", 
         paste(names(SummarizedExperiment::assays(se)), collapse = ", "))
  }
  
  beta_matrix <- SummarizedExperiment::assay(se, assay_name)
  
  # Convert to data frame format expected by EnsembleAge
  methylation_data <- data.frame(
    CGid = rownames(beta_matrix),
    beta_matrix,
    check.names = FALSE
  )
  
  # Handle sample sheet
  if (is.null(sample_sheet)) {
    if (verbose) cat("Extracting sample information from colData...\n")
    
    col_data <- as.data.frame(SummarizedExperiment::colData(se))
    
    # Ensure required columns exist
    if (!"Basename" %in% colnames(col_data)) {
      col_data$Basename <- colnames(se)
    }
    
    if (!"Age" %in% colnames(col_data)) {
      warning("Age column not found in colData. Using placeholder values.")
      col_data$Age <- rep(1.0, ncol(se))
    }
    
    sample_sheet <- col_data
  }
  
  if (verbose) {
    cat("Sample sheet columns:", paste(colnames(sample_sheet), collapse = ", "), "\n")
  }
  
  # Run age prediction using EnsembleAge
  results <- predict_all_clocks(methylation_data, sample_sheet, verbose = verbose)
  
  return(results)
}

#' Convert EnsembleAge Results to SummarizedExperiment
#'
#' Convert age prediction results back to a SummarizedExperiment object
#' for integration with other Bioconductor workflows.
#'
#' @param results Data frame from EnsembleAge prediction functions
#' @param original_se Optional original SummarizedExperiment to preserve metadata
#' @return SummarizedExperiment object with age predictions as rowData
#' @export
#' @examples
#' \donttest{
#' # Using simulated data
#' library(SummarizedExperiment)
#' 
#' # Create sample data
#' methylation_data <- data.frame(
#'   CGid = paste0("cg", 1:50),
#'   matrix(rbeta(50 * 4, 2, 2), nrow = 50, ncol = 4)
#' )
#' colnames(methylation_data)[-1] <- paste0("Sample_", 1:4)
#' 
#' sample_sheet <- data.frame(
#'   Basename = paste0("Sample_", 1:4),
#'   Age = c(0.5, 1.2, 2.1, 0.8)
#' )
#' 
#' # Get predictions
#' results <- predict_ensemble_static(methylation_data, sample_sheet, verbose = FALSE)
#' 
#' # Convert to SummarizedExperiment
#' se_results <- results_to_se(results)
#' }
results_to_se <- function(results, original_se = NULL) {
  
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is required but not installed.")
  }
  
  # Reshape results to wide format for SummarizedExperiment
  wide_results <- results %>%
    dplyr::select(Basename, epiClock, epiAge) %>%
    tidyr::pivot_wider(names_from = epiClock, values_from = epiAge, values_fn = mean)
  
  # Create assay matrix
  assay_matrix <- as.matrix(wide_results[, -1])
  rownames(assay_matrix) <- wide_results$Basename
  
  # Create column data with sample information
  # Select only columns that exist
  available_cols <- intersect(c("Basename", "Age", "Female", "Tissue", "SpeciesLatinName"), 
                             colnames(results))
  col_data <- results %>%
    dplyr::select(all_of(available_cols)) %>%
    dplyr::distinct()
  rownames(col_data) <- col_data$Basename
  
  # Create SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(predicted_ages = t(assay_matrix)),
    colData = col_data
  )
  
  return(se)
}

#' Extract Methylation Data from Common Bioconductor Objects
#'
#' Helper function to extract methylation beta values from various
#' Bioconductor objects commonly produced by sesame, minfi, etc.
#'
#' @param object Input object (SummarizedExperiment, matrix, or data.frame)
#' @param assay_name For SummarizedExperiment, name of assay to extract
#' @return Data frame in EnsembleAge format
#' @export
#' @examples
#' \donttest{
#' # With matrix
#' mat <- matrix(rbeta(100 * 4, 2, 2), nrow = 100, ncol = 4)
#' rownames(mat) <- paste0("cg", 1:100)
#' colnames(mat) <- paste0("Sample_", 1:4)
#' 
#' methylation_df <- extract_methylation_data(mat)
#' head(methylation_df)
#' }
extract_methylation_data <- function(object, assay_name = "beta") {
  
  if (methods::is(object, "SummarizedExperiment")) {
    # Extract from SummarizedExperiment
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment package is required")
    }
    
    beta_matrix <- SummarizedExperiment::assay(object, assay_name)
    methylation_data <- data.frame(
      CGid = rownames(beta_matrix),
      beta_matrix,
      check.names = FALSE
    )
    
  } else if (is.matrix(object)) {
    # Convert matrix to data frame
    methylation_data <- data.frame(
      CGid = rownames(object),
      object,
      check.names = FALSE
    )
    
  } else if (is.data.frame(object)) {
    # Assume already in correct format or convert
    if (!"CGid" %in% colnames(object) && !is.null(rownames(object))) {
      methylation_data <- data.frame(
        CGid = rownames(object),
        object,
        check.names = FALSE
      )
    } else {
      methylation_data <- object
    }
    
  } else {
    stop("Unsupported object type. Please provide SummarizedExperiment, matrix, or data.frame")
  }
  
  return(methylation_data)
}
