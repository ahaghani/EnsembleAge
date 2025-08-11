#' Create example methylation data
#'
#' This function creates example methylation data for testing and demonstration purposes.
#' The data mimics the structure of real methylation arrays with CpG sites and samples.
#' 
#' @param n_probes Number of CpG probes to include (default: 1000)
#' @param n_samples Number of samples to include (default: 20)
#' @param platform Platform type to simulate: "Mammal40k", "Mammal320k", or "Human450k" (default: "Mammal40k")
#' @param add_noise Logical indicating whether to add realistic noise (default: TRUE)
#' @return Data frame with CGid column and sample columns containing beta values
#' @export
#' @examples
#' # Create example data for Mammal40k platform
#' example_data <- create_example_methylation_data(n_probes = 500, n_samples = 10)
#' 
#' # Create larger dataset for Mammal320k
#' large_data <- create_example_methylation_data(n_probes = 5000, n_samples = 50, 
#'                                              platform = "Mammal320k")
create_example_methylation_data <- function(n_probes = 1000, n_samples = 20, 
                                           platform = "Mammal40k", add_noise = TRUE) {
  
  # Generate realistic CpG probe names based on platform
  if (platform == "Mammal40k" || platform == "Mammal320k") {
    cgids <- paste0("cg", sprintf("%08d", sample(10000000:99999999, n_probes, replace = FALSE)))
  } else if (platform == "Human450k" || platform == "HumanEPIC") {
    # Mix of cg and ch probes for human platforms
    n_cg <- round(n_probes * 0.95)
    n_ch <- n_probes - n_cg
    cg_probes <- paste0("cg", sprintf("%08d", sample(10000000:99999999, n_cg, replace = FALSE)))
    ch_probes <- paste0("ch", sprintf("%08d", sample(10000000:99999999, n_ch, replace = FALSE)))
    cgids <- c(cg_probes, ch_probes)
  } else {
    cgids <- paste0("cg", sprintf("%08d", sample(10000000:99999999, n_probes, replace = FALSE)))
  }
  
  # Generate sample names
  sample_names <- paste0("sample_", sprintf("%03d", 1:n_samples))
  
  # Create realistic methylation data
  # Use different distributions for different probe types
  methylation_data <- matrix(nrow = n_probes, ncol = n_samples)
  
  for (i in 1:n_probes) {
    # Different methylation patterns
    if (i <= n_probes * 0.3) {
      # Low methylation (promoter regions)
      base_values <- rbeta(n_samples, 1, 9)  # Biased towards low values
    } else if (i <= n_probes * 0.6) {
      # Medium methylation
      base_values <- rbeta(n_samples, 2, 2)  # More uniform
    } else {
      # High methylation (gene bodies, intergenic)
      base_values <- rbeta(n_samples, 9, 1)  # Biased towards high values
    }
    
    # Add age-related pattern to some probes
    if (i <= n_probes * 0.1) {
      age_effect <- seq(-0.2, 0.2, length.out = n_samples)
      base_values <- pmin(pmax(base_values + age_effect, 0), 1)
    }
    
    methylation_data[i, ] <- base_values
  }
  
  # Add noise if requested
  if (add_noise) {
    noise <- matrix(rnorm(n_probes * n_samples, 0, 0.05), nrow = n_probes)
    methylation_data <- pmin(pmax(methylation_data + noise, 0), 1)
  }
  
  # Create data frame
  result <- data.frame(
    CGid = cgids,
    methylation_data,
    stringsAsFactors = FALSE
  )
  names(result)[-1] <- sample_names
  
  return(result)
}

#' Create example sample information
#'
#' This function creates example sample information data that matches the structure
#' expected by EnsembleAge prediction functions.
#' 
#' @param sample_names Character vector of sample names (default: NULL, will generate)
#' @param n_samples Number of samples if sample_names not provided (default: 20)
#' @param species Species for samples (default: "Mus musculus")
#' @param add_groups Logical indicating whether to add experimental groups (default: TRUE)
#' @return Data frame with sample information
#' @export
#' @examples
#' # Create sample info for 10 samples
#' sample_info <- create_example_sample_info(n_samples = 10)
#' 
#' # Create sample info for specific sample names
#' my_samples <- paste0("sample_", 1:5)
#' sample_info <- create_example_sample_info(sample_names = my_samples, 
#'                                          species = "Rattus norvegicus")
create_example_sample_info <- function(sample_names = NULL, n_samples = 20, 
                                      species = "Mus musculus", add_groups = TRUE) {
  
  if (is.null(sample_names)) {
    sample_names <- paste0("sample_", sprintf("%03d", 1:n_samples))
  } else {
    n_samples <- length(sample_names)
  }
  
  # Generate realistic ages based on species
  if (species == "Mus musculus") {
    # Mouse ages in years (0.1 to 3 years)
    ages <- round(runif(n_samples, 0.1, 3), 2)
  } else if (species == "Rattus norvegicus") {
    # Rat ages in years (0.1 to 3.5 years)
    ages <- round(runif(n_samples, 0.1, 3.5), 2)
  } else if (species == "Homo sapiens") {
    # Human ages in years (1 to 90 years)
    ages <- round(runif(n_samples, 1, 90), 0)
  } else {
    # Default to mouse-like ages
    ages <- round(runif(n_samples, 0.1, 3), 2)
  }
  
  # Generate sex information
  female_status <- sample(c(0, 1), n_samples, replace = TRUE)
  
  # Generate tissue types
  if (species %in% c("Mus musculus", "Rattus norvegicus")) {
    tissues <- sample(c("Blood", "Liver", "Heart", "Kidney", "Muscle", "Brain", "Cortex", "Tail"), 
                     n_samples, replace = TRUE)
  } else {
    tissues <- sample(c("Blood", "Saliva", "Buccal", "Brain", "Liver"), 
                     n_samples, replace = TRUE)
  }
  
  # Create basic sample info
  sample_info <- data.frame(
    Basename = sample_names,
    Age = ages,
    SpeciesLatinName = species,
    Female = female_status,
    Tissue = tissues,
    stringsAsFactors = FALSE
  )
  
  # Add experimental groups if requested
  if (add_groups) {
    n_groups <- min(4, max(2, n_samples %/% 5))
    sample_info$Group <- sample(paste0("Group_", LETTERS[1:n_groups]), n_samples, replace = TRUE)
    
    # Add treatment information
    sample_info$Treatment <- sample(c("Control", "Treatment"), n_samples, replace = TRUE)
  }
  
  return(sample_info)
}

#' Run complete example analysis
#'
#' This function demonstrates a complete analysis workflow using example data.
#' It creates synthetic data, runs preprocessing, and performs age predictions.
#' 
#' @param n_probes Number of probes in example data (default: 2000)
#' @param n_samples Number of samples in example data (default: 30)
#' @param platform Platform to simulate (default: "Mammal40k")
#' @param verbose Logical indicating whether to print progress (default: TRUE)
#' @return List containing example data and results
#' @export
#' @examples
#' \dontrun{
#' # Run a complete example analysis
#' example_results <- run_example_analysis(n_probes = 1000, n_samples = 20)
#' 
#' # View the results
#' head(example_results$age_predictions)
#' 
#' # Check data quality
#' example_results$data_validation
#' }
run_example_analysis <- function(n_probes = 2000, n_samples = 30, 
                                platform = "Mammal40k", verbose = TRUE) {
  
  if (verbose) cat("Creating example methylation data...\n")
  
  # Create example data
  methylation_data <- create_example_methylation_data(
    n_probes = n_probes, 
    n_samples = n_samples, 
    platform = platform
  )
  
  sample_info <- create_example_sample_info(
    sample_names = names(methylation_data)[-1]
  )
  
  if (verbose) {
    cat("Created data with", nrow(methylation_data), "probes and", 
        nrow(sample_info), "samples\n")
  }
  
  # Validate data
  if (verbose) cat("Validating data format...\n")
  validation <- validate_data_format(methylation_data, sample_info)
  
  if (verbose) {
    cat("Platform detected:", validation$platform, "\n")
    cat("Data validation:", ifelse(validation$valid, "PASSED", "FAILED"), "\n")
  }
  
  # Check probe coverage
  if (verbose) cat("Checking probe coverage for clocks...\n")
  coverage <- check_probe_coverage(methylation_data)
  
  if (verbose) {
    cat("Top 5 clocks by probe coverage:\n")
    print(head(coverage, 5))
  }
  
  # Run preprocessing
  if (verbose) cat("Preprocessing data...\n")
  processed_data <- preprocess_methylation_data(
    methylation_data, 
    sample_info, 
    verbose = verbose
  )
  
  # Run age prediction
  if (verbose) cat("Running age predictions...\n")
  age_predictions <- predictAgeAndAgeAcc(
    processed_data$data, 
    processed_data$samples
  )
  
  if (verbose) {
    cat("Prediction complete! Predicted", ncol(age_predictions) - 1, "age-related measures\n")
    cat("Available predictions:", paste(names(age_predictions)[grepl("epiAge", names(age_predictions))][1:3], collapse = ", "), "...\n")
  }
  
  return(list(
    methylation_data = methylation_data,
    sample_info = sample_info,
    processed_data = processed_data,
    data_validation = validation,
    probe_coverage = coverage,
    age_predictions = age_predictions
  ))
}
