#' Detect and organize FASTQ files
#'
#' This function scans a directory for FASTQ files and organizes them into a structured
#' format for downstream analysis. It can detect paired-end files based on naming patterns.
#'
#' @param input_dir Character string specifying the directory containing FASTQ files
#' @param pattern Character string specifying a pattern to match in file names (default: ".fastq|.fq")
#' @param paired Logical indicating whether to look for paired-end files (default: TRUE)
#' @param forward_pattern Character string pattern to identify forward reads (default: "_R1")
#' @param reverse_pattern Character string pattern to identify reverse reads (default: "_R2")
#' @param sample_pattern Regular expression to extract sample names (default: "(.+?)(_S[0-9]+)?(_L[0-9]+)?(_R[12])")
#' @param recursive Logical indicating whether to search directories recursively (default: FALSE)
#'
#' @return A tibble with columns for sample names, forward read paths, reverse read paths (if paired), 
#'         and file sizes
#'
#' @export
#' @importFrom dplyr filter mutate select arrange left_join %>%
#' @importFrom tibble tibble
#' @importFrom readr read_lines
#' @importFrom stringr str_extract str_replace_all str_detect
#' @importFrom purrr map_dbl
#'
#' @examples
#' \dontrun{
#' # For paired-end files
#' fastq_files <- detectFastqFiles("path/to/fastq/dir")
#' 
#' # For single-end files
#' fastq_files <- detectFastqFiles("path/to/fastq/dir", paired = FALSE)
#' }
detectFastqFiles <- function(input_dir, 
                            pattern = ".fastq|.fq",
                            paired = TRUE,
                            forward_pattern = "_R1",
                            reverse_pattern = "_R2",
                            sample_pattern = "(.+?)(_S[0-9]+)?(_L[0-9]+)?(_R[12])",
                            recursive = FALSE) {
  
  # Input validation
  if (!._dirExists(input_dir)) {
    stop("Input directory does not exist: ", input_dir)
  }
  
  # Find all files matching pattern
  all_files <- list.files(path = input_dir, 
                         pattern = pattern, 
                         full.names = TRUE, 
                         recursive = recursive)
  
  if (length(all_files) == 0) {
    stop("No files matching pattern found in: ", input_dir)
  }
  
  # Extract file information
  file_info <- tibble(
    file_path = all_files,
    file_name = basename(file_path),
    file_size = file.size(file_path),
    sample_name = str_extract(file_name, sample_pattern, group = 1)
  )
  
  if (any(is.na(file_info$sample_name))) {
    warning("Could not extract sample names for some files using pattern: ", sample_pattern)
  }
  
  # Handle paired or single-end files
  if (paired) {
    forward_files <- file_info %>%
      filter(str_detect(file_name, forward_pattern)) %>%
      select(sample_name, forward = file_path, forward_size = file_size)
    
    reverse_files <- file_info %>%
      filter(str_detect(file_name, reverse_pattern)) %>%
      select(sample_name, reverse = file_path, reverse_size = file_size)
    
    # Join forward and reverse files by sample name
    fastq_data <- forward_files %>%
      left_join(reverse_files, by = "sample_name") %>%
      arrange(sample_name)
    
    # Check if all forward files have matching reverse files
    missing_reverse <- is.na(fastq_data$reverse)
    if (any(missing_reverse)) {
      warning("Missing reverse files for samples: ", 
              paste(fastq_data$sample_name[missing_reverse], collapse = ", "))
    }
  } else {
    # Single-end mode
    fastq_data <- file_info %>%
      select(sample_name, file_path, file_size) %>%
      arrange(sample_name)
  }
  
  # Provide summary
  message("Detected ", nrow(fastq_data), " samples")
  
  return(fastq_data)
}

#' Validate FASTQ files
#'
#' This function performs basic validation on FASTQ files to ensure they are readable
#' and contain valid FASTQ format.
#'
#' @param fastq_files Tibble with FASTQ file information as produced by detectFastqFiles()
#' @param check_lines Number of lines to check from beginning of each file (default: 40)
#' @param paired Logical indicating whether files are paired (default: TRUE)
#'
#' @return The input tibble with an additional column indicating validation status
#'
#' @export
#' @importFrom dplyr mutate
#' @importFrom purrr map_lgl
#' @importFrom readr read_lines
#'
#' @examples
#' \dontrun{
#' fastq_files <- detectFastqFiles("path/to/fastq/dir")
#' validated_files <- validateFastqFiles(fastq_files)
#' }
validateFastqFiles <- function(fastq_files, check_lines = 40, paired = TRUE) {
  
  # Function to check if a file is valid FASTQ
  ._isValidFastq <- function(file_path, check_lines = 40) {
    if (!._fileExists(file_path)) {
      return(FALSE)
    }
    
    tryCatch({
      # Read the first few lines
      lines <- read_lines(file_path, n_max = check_lines)
      
      if (length(lines) < 4) {
        return(FALSE)
      }
      
      # Check FASTQ format: header (@), sequence, separator (+), quality
      for (i in seq(1, min(length(lines), check_lines), by = 4)) {
        if (i + 3 > length(lines)) break
        
        if (!str_detect(lines[i], "^@") || 
            !str_detect(lines[i+2], "^\\+")) {
          return(FALSE)
        }
        
        # Check that sequence and quality lines have same length
        if (nchar(lines[i+1]) != nchar(lines[i+3])) {
          return(FALSE)
        }
      }
      
      return(TRUE)
    }, error = function(e) {
      message("Error validating file: ", file_path, " - ", e$message)
      return(FALSE)
    })
  }
  
  # Validate forward files
  fastq_files <- fastq_files %>%
    mutate(forward_valid = map_lgl(forward, ._isValidFastq, check_lines = check_lines))
  
  # Validate reverse files if paired
  if (paired && "reverse" %in% names(fastq_files)) {
    fastq_files <- fastq_files %>%
      mutate(reverse_valid = map_lgl(reverse, ._isValidFastq, check_lines = check_lines))
  }
  
  # Report validation results
  invalid_count <- sum(!fastq_files$forward_valid)
  if (paired && "reverse_valid" %in% names(fastq_files)) {
    invalid_count <- invalid_count + sum(!fastq_files$reverse_valid, na.rm = TRUE)
  }
  
  if (invalid_count > 0) {
    warning("Found ", invalid_count, " invalid FASTQ files")
  } else {
    message("All FASTQ files validated successfully")
  }
  
  return(fastq_files)
}