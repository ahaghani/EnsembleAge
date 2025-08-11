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
prepare_sample_sheet <- function(samps, verbose = TRUE) {
  
  missing_vars <- character(0)
  defaults_applied <- character(0)
  
  # Ensure Basename column exists
  if (!"Basename" %in% names(samps)) {
    if ("Sample_Name" %in% names(samps)) {
      samps$Basename <- samps$Sample_Name
      if (verbose) cat("Using 'Sample_Name' column as 'Basename'\n")
    } else {
      stop("Sample sheet must contain either 'Basename' or 'Sample_Name' column")
    }
  }
  
  # Ensure Age column exists - this is truly required for clock predictions
  if (!"Age" %in% names(samps)) {
    stop("Sample sheet must contain 'Age' column with chronological ages in years")
  }
  
  # Add missing optional columns with defaults and track what was added
  if (!"Female" %in% names(samps)) {
    samps <- samps %>% mutate(Female = NA)
    missing_vars <- c(missing_vars, "Female")
    defaults_applied <- c(defaults_applied, "Female = NA (sex unknown)")
  }
  
  if (!"Tissue" %in% names(samps)) {
    samps <- samps %>% mutate(Tissue = "Unknown")
    missing_vars <- c(missing_vars, "Tissue")
    defaults_applied <- c(defaults_applied, "Tissue = 'Unknown'")
  }
  
  if (!"SpeciesLatinName" %in% names(samps)) {
    samps <- samps %>% mutate(SpeciesLatinName = "Mus musculus")
    missing_vars <- c(missing_vars, "SpeciesLatinName")
    defaults_applied <- c(defaults_applied, "SpeciesLatinName = 'Mus musculus'")
  }
  
  # Notify user about defaults applied
  if (length(missing_vars) > 0 && verbose) {
    cat("Missing variables detected. Applied defaults:\n")
    for (default in defaults_applied) {
      cat("  -", default, "\n")
    }
    cat("For better predictions, consider providing these variables in your sample sheet.\n")
  }
  
  # Convert data types with error handling
  tryCatch({
    samps$Age <- as.numeric(samps$Age)
    if (any(is.na(samps$Age))) {
      warning("Some Age values could not be converted to numeric. Check your data.")
    }
  }, error = function(e) {
    stop("Error converting Age column to numeric: ", e$message)
  })
  
  # Handle Female column conversion
  if (!all(is.na(samps$Female))) {
    tryCatch({
      samps$Female <- as.numeric(samps$Female)
      # Validate Female values (should be 0, 1, or NA)
      invalid_female <- samps$Female[!is.na(samps$Female) & !samps$Female %in% c(0, 1)]
      if (length(invalid_female) > 0) {
        warning("Female column should contain only 0 (male), 1 (female), or NA. Invalid values found.")
      }
    }, error = function(e) {
      warning("Error converting Female column: ", e$message, ". Setting to NA.")
      samps$Female <- NA
    })
  }
  
  # Ensure character columns
  samps$Tissue <- as.character(samps$Tissue)
  samps$SpeciesLatinName <- as.character(samps$SpeciesLatinName)
  samps$Basename <- as.character(samps$Basename)
  
  # Validate that we have samples
  if (nrow(samps) == 0) {
    stop("Sample sheet is empty")
  }
  
  # Check for duplicate Basenames
  if (any(duplicated(samps$Basename))) {
    duplicates <- samps$Basename[duplicated(samps$Basename)]
    warning("Duplicate sample names found: ", paste(unique(duplicates), collapse = ", "))
  }
  
  if (verbose) {
    cat("Sample sheet prepared:", nrow(samps), "samples\n")
  }
  
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

#' Map probe IDs to CpG IDs for Mammal320k data
#'
#' This function maps Mammal320k probe IDs (with suffixes like _BC21) to standard CpG IDs
#' using the annotation file.
#' 
#' @param dat0sesame Data frame containing methylation data with probe IDs
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @return Data frame with CpG IDs mapped and duplicates handled
#' @export
#' @examples
#' \dontrun{
#' mapped_data <- map_mammal320k_probes(methylation_data)
#' }
map_mammal320k_probes <- function(dat0sesame, verbose = TRUE) {
  
  # Load Mammal320k annotation
  annotation_file <- system.file("data", "Mus musculus. Mammalian 320k. mm10.Amin.V10.RDS", 
                                 package = "EnsembleAge")
  if (annotation_file == "") {
    annotation_file <- file.path("data", "Mus musculus. Mammalian 320k. mm10.Amin.V10.RDS")
  }
  
  if (!file.exists(annotation_file)) {
    warning("Mammal320k annotation file not found. Cannot map probe IDs.")
    return(dat0sesame)
  }
  
  annotation <- readRDS(annotation_file)
  
  if (verbose) {
    cat("Loaded Mammal320k annotation with", nrow(annotation), "probes\n")
  }
  
  # Check if CGid column exists, if not assume first column contains probe IDs
  if ("CGid" %in% names(dat0sesame)) {
    probe_col <- "CGid"
  } else {
    probe_col <- names(dat0sesame)[1]
    if (verbose) cat("No CGid column found, using", probe_col, "as probe ID column\n")
  }
  
  # Get current probe IDs
  current_probes <- dat0sesame[[probe_col]]
  
  # Check if probes need mapping (have suffixes)
  has_suffixes <- sum(grepl("_[A-Z]+[0-9]*$", current_probes))
  
  if (has_suffixes == 0) {
    if (verbose) cat("Probe IDs appear to be clean CpG IDs already. No mapping needed.\n")
    if (probe_col != "CGid") {
      names(dat0sesame)[names(dat0sesame) == probe_col] <- "CGid"
    }
    return(dat0sesame)
  }
  
  if (verbose) {
    cat("Found", has_suffixes, "probes with suffixes that need mapping\n")
  }
  
  # Create mapping from annotation
  probe_mapping <- annotation %>%
    dplyr::select(Probe_ID, CGid) %>%
    filter(!is.na(CGid) & !is.na(Probe_ID))
  
  # Map probe IDs to CpG IDs
  dat_mapped <- dat0sesame %>%
    dplyr::rename(Probe_ID = !!probe_col) %>%
    left_join(probe_mapping, by = "Probe_ID") %>%
    filter(!is.na(CGid)) %>%
    dplyr::select(-Probe_ID)
  
  # Handle duplicates by taking the mean
  if (any(duplicated(dat_mapped$CGid))) {
    if (verbose) cat("Found duplicate CpG IDs after mapping. Taking mean values...\n")
    
    sample_cols <- names(dat_mapped)[names(dat_mapped) != "CGid"]
    dat_mapped <- dat_mapped %>%
      group_by(CGid) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(sample_cols), mean, na.rm = TRUE), .groups = "drop")
  }
  
  if (verbose) {
    cat("Mapping complete:", nrow(dat_mapped), "unique CpG IDs retained\n")
    cat("Original probes:", nrow(dat0sesame), "-> Mapped CpGs:", nrow(dat_mapped), "\n")
  }
  
  return(dat_mapped)
}
