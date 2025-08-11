#' Detect methylation array platform
#'
#' This function detects the methylation array platform based on the number of CpG sites
#' and specific probe patterns in the data.
#' 
#' @param dat0sesame Data frame containing methylation data with CGid column
#' @return Character string indicating the detected platform: "Mammal320k", "Mammal40k", "Human450k", "HumanEPIC", or "Unknown"
#' @export
#' @examples
#' \dontrun{
#' # Example with methylation data
#' platform <- detect_platform(methylation_data)
#' }
detect_platform <- function(dat0sesame) {
  if (!"CGid" %in% names(dat0sesame)) {
    stop("Data must contain a 'CGid' column with CpG identifiers")
  }
  
  n_probes <- nrow(dat0sesame)
  probe_ids <- dat0sesame$CGid
  
  # Check for specific platform patterns
  if (n_probes > 250000) {
    return("Mammal320k")
  } else if (n_probes > 30000 && n_probes < 50000) {
    return("Mammal40k")
  } else if (n_probes > 400000 && n_probes < 500000) {
    return("Human450k")
  } else if (n_probes > 800000) {
    return("HumanEPIC")
  } else {
    # Additional checks based on probe naming patterns
    if (any(grepl("^cg", probe_ids)) && any(grepl("^ch", probe_ids))) {
      if (n_probes > 800000) {
        return("HumanEPIC")
      } else {
        return("Human450k")
      }
    } else if (any(grepl("^cg", probe_ids))) {
      return("Mammal40k")
    } else {
      return("Unknown")
    }
  }
}

#' Validate methylation data format
#'
#' This function validates that the methylation data is in the correct format
#' for EnsembleAge predictions.
#' 
#' @param dat0sesame Data frame containing methylation data
#' @param samps Data frame containing sample information
#' @return List with validation results and suggestions
#' @export
#' @examples
#' \dontrun{
#' validation <- validate_data_format(methylation_data, sample_info)
#' if (!validation$valid) {
#'   stop(validation$message)
#' }
#' }
validate_data_format <- function(dat0sesame, samps) {
  issues <- character(0)
  suggestions <- character(0)
  
  # Check dat0sesame format
  if (!"CGid" %in% names(dat0sesame)) {
    issues <- c(issues, "dat0sesame must contain a 'CGid' column")
    suggestions <- c(suggestions, "Add a 'CGid' column with CpG identifiers")
  }
  
  if (ncol(dat0sesame) < 2) {
    issues <- c(issues, "dat0sesame must contain sample columns in addition to CGid")
  }
  
  # Check for numeric methylation values
  numeric_cols <- sapply(dat0sesame[, !names(dat0sesame) %in% "CGid"], is.numeric)
  if (!all(numeric_cols)) {
    issues <- c(issues, "All sample columns must contain numeric methylation values")
  }
  
  # Check value ranges
  if (any(sapply(dat0sesame[, !names(dat0sesame) %in% "CGid"], function(x) any(x < 0 | x > 1, na.rm = TRUE)))) {
    issues <- c(issues, "Methylation values should be between 0 and 1 (beta values)")
    suggestions <- c(suggestions, "Convert M-values to beta values if necessary")
  }
  
  # Check samps format
  if (!"Basename" %in% names(samps)) {
    issues <- c(issues, "samps must contain a 'Basename' column with sample identifiers")
  }
  
  # Check if sample names match
  sample_cols <- names(dat0sesame)[!names(dat0sesame) %in% "CGid"]
  if ("Basename" %in% names(samps)) {
    missing_samples <- setdiff(samps$Basename, sample_cols)
    if (length(missing_samples) > 0) {
      issues <- c(issues, paste("Some samples in samps$Basename are not found in dat0sesame:", 
                                paste(missing_samples[1:min(5, length(missing_samples))], collapse = ", ")))
    }
  }
  
  # Platform detection
  platform <- detect_platform(dat0sesame)
  
  return(list(
    valid = length(issues) == 0,
    platform = platform,
    issues = issues,
    suggestions = suggestions,
    n_probes = nrow(dat0sesame),
    n_samples = ncol(dat0sesame) - 1
  ))
}

#' Prepare sample sheet with default values
#'
#' This function prepares the sample sheet by adding default values for missing columns
#' required by the EnsembleAge prediction functions.
#' 
#' @param samps Data frame containing sample information
#' @param default_species Character string for default species (default: "Mus musculus")
#' @param default_age Numeric value for default age (default: 0)
#' @return Data frame with standardized sample information
#' @export
#' @importFrom dplyr mutate
#' @examples
#' \dontrun{
#' sample_info <- data.frame(Basename = c("sample1", "sample2"))
#' prepared_samples <- prepare_sample_sheet(sample_info)
#' }
prepare_sample_sheet <- function(samps, default_species = "Mus musculus", default_age = 0) {
  
  # Add default columns if missing
  if (!"SpeciesLatinName" %in% names(samps)) {
    samps <- samps %>% mutate(SpeciesLatinName = default_species)
  }
  
  if (!"Age" %in% names(samps)) {
    samps <- samps %>% mutate(Age = default_age)
  }
  
  if (!"Female" %in% names(samps)) {
    samps <- samps %>% mutate(Female = NA)
  }
  
  if (!"Tissue" %in% names(samps)) {
    samps <- samps %>% mutate(Tissue = NA)
  }
  
  # Clean up existing values
  samps <- samps %>% 
    mutate(Age = ifelse(is.na(Age), default_age, Age)) %>% 
    mutate(SpeciesLatinName = ifelse(is.na(SpeciesLatinName), default_species, SpeciesLatinName))
  
  return(samps)
}

#' Get available clocks for a platform
#'
#' This function returns the available clocks for a specific methylation platform.
#' 
#' @param platform Character string indicating the platform
#' @return Character vector of available clock names
#' @export
#' @examples
#' available_clocks <- get_available_clocks("Mammal40k")
get_available_clocks <- function(platform) {
  clock_file_path <- system.file("data", "Clock coefficients.RDS", package = "EnsembleAge")
  if (clock_file_path == "") {
    clock_file_path <- file.path("data", "Clock coefficients.RDS")
  }
  
  if (!file.exists(clock_file_path)) {
    warning("Clock coefficients file not found. Please ensure the package is properly installed.")
    return(character(0))
  }
  
  epiclocks <- readRDS(clock_file_path)
  
  # Return all available clock families and individual clocks
  all_clocks <- character(0)
  for (family_name in names(epiclocks)) {
    family_clocks <- epiclocks[[family_name]]
    clock_names <- paste(family_name, names(family_clocks), sep = ".")
    all_clocks <- c(all_clocks, clock_names)
  }
  
  return(all_clocks)
}

#' Check probe coverage for platform
#'
#' This function checks what percentage of required probes are available in the input data
#' for each clock type.
#' 
#' @param dat0sesame Data frame containing methylation data with CGid column
#' @return Data frame with clock names and probe coverage percentages
#' @export
#' @examples
#' \dontrun{
#' coverage <- check_probe_coverage(methylation_data)
#' }
check_probe_coverage <- function(dat0sesame) {
  if (!"CGid" %in% names(dat0sesame)) {
    stop("Data must contain a 'CGid' column")
  }
  
  clock_file_path <- system.file("data", "Clock coefficients.RDS", package = "EnsembleAge")
  if (clock_file_path == "") {
    clock_file_path <- file.path("data", "Clock coefficients.RDS")
  }
  
  if (!file.exists(clock_file_path)) {
    warning("Clock coefficients file not found.")
    return(data.frame(Clock = character(0), Coverage = numeric(0)))
  }
  
  epiclocks <- readRDS(clock_file_path)
  available_probes <- dat0sesame$CGid
  
  coverage_results <- data.frame(
    Clock = character(0),
    Required_Probes = integer(0),
    Available_Probes = integer(0),
    Coverage_Percent = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (family_name in names(epiclocks)) {
    family_clocks <- epiclocks[[family_name]]
    for (clock_name in names(family_clocks)) {
      clock <- family_clocks[[clock_name]]
      required_probes <- clock$CGid
      available_in_data <- sum(required_probes %in% available_probes)
      coverage_percent <- round(available_in_data / length(required_probes) * 100, 1)
      
      coverage_results <- rbind(coverage_results, data.frame(
        Clock = paste(family_name, clock_name, sep = "."),
        Required_Probes = length(required_probes),
        Available_Probes = available_in_data,
        Coverage_Percent = coverage_percent,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(coverage_results[order(coverage_results$Coverage_Percent, decreasing = TRUE), ])
}
