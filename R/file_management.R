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
  
  # Return detected files (no message output)
  
  return(fastq_data)
}

#' Detect and organize FASTQ files with advanced matching capabilities
#'
#' This function provides enhanced file detection with fuzzy matching, metadata extraction,
#' Bayesian confidence scoring, directory structure analysis, and content verification.
#'
#' @param input_dir Character string specifying the directory containing FASTQ files
#' @param pattern Character string specifying a pattern to match in file names (default: ".fastq|.fq")
#' @param paired Logical indicating whether to look for paired-end files (default: TRUE)
#' @param fuzzy_threshold Numeric value between 0 and 1 indicating the maximum normalized
#'        string distance to consider files as potential matches (default: 0.2)
#' @param extract_metadata Logical indicating whether to attempt to extract metadata from
#'        filenames (default: TRUE)
#' @param verify_content Logical indicating whether to verify file content (default: TRUE)
#' @param recursive Logical indicating whether to search directories recursively (default: TRUE)
#' @param confidence_threshold Numeric value between 0 and 1 indicating the minimum confidence
#'        score to include a match in the results (default: 0.7)
#' @param platform Character string specifying the sequencing platform, or "auto" to attempt
#'        to detect (options: "auto", "illumina", "pacbio", "nanopore", "iontorrent") (default: "auto")
#' @param verify_lines Integer specifying the number of lines to check when verifying content (default: 1000)
#' @param known_patterns List containing named elements with platform-specific naming patterns.
#'        If NULL, uses default patterns. (default: NULL)
#'
#' @return A tibble with columns for sample names, forward and reverse read paths, confidence
#'         scores, metadata, and other information
#'
#' @export
#' @importFrom dplyr filter mutate select arrange left_join group_by ungroup %>% 
#' @importFrom tibble tibble
#' @importFrom readr read_lines
#' @importFrom stringr str_extract str_extract_all str_replace_all str_detect str_match str_split
#' @importFrom purrr map map_dbl map_lgl map2 map2_dbl pmap
#' @importFrom tidyr unnest
#' @importFrom stats setNames
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' fastq_files <- detectFastqFilesAdvanced("path/to/fastq/dir")
#'
#' # With lower confidence threshold for ambiguous naming
#' fastq_files <- detectFastqFilesAdvanced("path/to/fastq/dir", confidence_threshold = 0.5)
#'
#' # For a specific platform with content verification
#' fastq_files <- detectFastqFilesAdvanced(
#'   "path/to/fastq/dir", 
#'   platform = "pacbio", 
#'   verify_content = TRUE
#' )
#' }
detectFastqFilesAdvanced <- function(input_dir,
                                    pattern = ".fastq|.fq",
                                    paired = TRUE,
                                    fuzzy_threshold = 0.2,
                                    extract_metadata = TRUE,
                                    verify_content = TRUE,
                                    recursive = TRUE,
                                    confidence_threshold = 0.7,
                                    platform = "auto",
                                    verify_lines = 1000,
                                    known_patterns = NULL) {
  
  # Input validation
  if (!._dirExists(input_dir)) {
    stop("Input directory does not exist: ", input_dir)
  }
  
  # Check for required packages
  if (!requireNamespace("stringdist", quietly = TRUE)) {
    stop("Package 'stringdist' is required for fuzzy matching. Please install it with install.packages('stringdist')")
  }
  
  # Get all FASTQ files from the directory structure
  message("Scanning for FASTQ files...")
  file_list <- ._findFastqFiles(input_dir, pattern, recursive)
  if (length(file_list) == 0) {
    stop("No FASTQ files found in the specified directory")
  }
  message("Found ", length(file_list), " FASTQ files")
  
  # Analyze directory structure for organizational patterns
  dir_structure <- ._analyzeDirectoryStructure(file_list)
  
  # Detect platform if set to auto
  if (platform == "auto") {
    platform <- ._detectSequencingPlatform(file_list, dir_structure)
    message("Detected platform: ", platform)
  }
  
  # Get platform-specific patterns
  patterns <- ._getPlatformPatterns(platform, known_patterns)
  
  # Extract file information and metadata
  file_info <- ._extractFileInfo(file_list, patterns, extract_metadata)
  
  # Pair files with confidence scores
  if (paired) {
    message("Pairing files using fuzzy matching...")
    paired_files <- ._pairFastqFiles(
      file_info, 
      patterns, 
      fuzzy_threshold, 
      platform
    )
    
    # Filter by confidence threshold
    filtered_pairs <- paired_files %>%
      filter(confidence >= confidence_threshold)
    
    if (nrow(filtered_pairs) == 0) {
      warning("No file pairs met the confidence threshold. Consider lowering the threshold.")
      # Return empty result with proper structure
      return(tibble(
        sample_name = character(),
        forward = character(),
        reverse = character(),
        confidence = numeric(),
        metadata = list(),
        platform = character()
      ))
    }
    
    result <- filtered_pairs
  } else {
    # Single-end mode
    result <- file_info %>%
      mutate(
        confidence = 1.0,
        platform = platform
      ) %>%
      select(sample_name, file_path, file_size, confidence, metadata, platform)
  }
  
  # Verify content if requested
  if (verify_content) {
    message("Verifying file content...")
    result <- ._verifyFastqContent(result, paired, verify_lines)
  }
  
  # Final organization of results
  if (paired) {
    organized_result <- result %>%
      arrange(desc(confidence), sample_name)
  } else {
    organized_result <- result %>%
      arrange(sample_name)
  }
  
  # Report results
  message("Found ", nrow(organized_result), " samples with ", 
          ifelse(paired, "paired", "single"), " reads")
  
  return(organized_result)
}

#' Find FASTQ files in a directory
#'
#' @param input_dir Character string specifying the directory to search
#' @param pattern Character string specifying the pattern to match FASTQ files
#' @param recursive Logical indicating whether to search recursively
#'
#' @return Character vector of file paths
#' @keywords internal
._findFastqFiles <- function(input_dir, pattern, recursive) {
  # Get all files
  files <- list.files(
    path = input_dir,
    pattern = pattern,
    full.names = TRUE,
    recursive = recursive
  )
  
  # Filter to keep only files (not directories)
  files[file.info(files)$isdir == FALSE]
}

#' Analyze directory structure for organizational patterns
#'
#' @param file_list Character vector of file paths
#'
#' @return List with directory structure information
#' @keywords internal
._analyzeDirectoryStructure <- function(file_list) {
  # Extract directories from file paths
  dirs <- unique(dirname(file_list))
  
  # Analyze directory depth
  depths <- sapply(strsplit(dirs, "/"), length)
  
  # Analyze directory naming patterns
  dir_names <- basename(dirs)
  
  # Look for common organizational patterns
  potential_batch_dirs <- grep("^(batch|run|plate|exp)", dir_names, value = TRUE, ignore.case = TRUE)
  potential_date_dirs <- grep("^[0-9]{4}[-_]?[0-9]{2}[-_]?[0-9]{2}", dir_names, value = TRUE)
  
  # Count files per directory
  files_per_dir <- table(dirname(file_list))
  
  # Analyze nested structure (for multi-timepoint, multi-platform studies)
  parent_dirs <- unique(dirname(dirs))
  nested_structure <- length(parent_dirs) < length(dirs)
  
  # Return structure analysis
  list(
    total_dirs = length(dirs),
    max_depth = max(depths),
    avg_depth = mean(depths),
    potential_batch_dirs = potential_batch_dirs,
    potential_date_dirs = potential_date_dirs,
    files_per_dir = files_per_dir,
    nested = nested_structure,
    parent_dirs = parent_dirs
  )
}

#' Detect sequencing platform from file naming conventions
#'
#' @param file_list Character vector of file paths
#' @param dir_structure List with directory structure information
#'
#' @return Character string indicating the detected platform
#' @keywords internal
#' @importFrom stringr str_detect
._detectSequencingPlatform <- function(file_list, dir_structure) {
  # Extract file names
  file_names <- basename(file_list)
  
  # Platform detection patterns
  illumina_pattern <- "_S[0-9]+_L[0-9]{3}_R[12]_001|_R[12]_001\\.fastq|\\.illumina\\."
  pacbio_pattern <- "\\.subreads\\.bam|\\.ccs\\.bam|\\.hifi_reads\\.fastq|m[0-9]{5}_[0-9]{6}"
  nanopore_pattern <- "barcode[0-9]{2}|NB[0-9]{2}|pass/|fast5|ONT"
  iontorrent_pattern <- "IonXpress|Ion_Torrent|PGMID"
  
  # Count matches for each platform
  illumina_count <- sum(str_detect(file_names, illumina_pattern))
  pacbio_count <- sum(str_detect(file_names, pacbio_pattern))
  nanopore_count <- sum(str_detect(file_names, nanopore_pattern))
  iontorrent_count <- sum(str_detect(file_names, iontorrent_pattern))
  
  # Also check directory names for platform hints
  dir_names <- c(basename(dir_structure$parent_dirs), basename(dirname(file_list)))
  
  illumina_dir_count <- sum(str_detect(dir_names, "illumina|miseq|nextseq|hiseq|novaseq", ignore.case = TRUE))
  pacbio_dir_count <- sum(str_detect(dir_names, "pacbio|sequel|rsii", ignore.case = TRUE))
  nanopore_dir_count <- sum(str_detect(dir_names, "nanopore|minion|gridion|promethion", ignore.case = TRUE))
  iontorrent_dir_count <- sum(str_detect(dir_names, "iontorrent|ion_torrent|pgm|s5", ignore.case = TRUE))
  
  # Combine file and directory evidence with weights
  illumina_evidence <- illumina_count * 2 + illumina_dir_count * 3
  pacbio_evidence <- pacbio_count * 2 + pacbio_dir_count * 3
  nanopore_evidence <- nanopore_count * 2 + nanopore_dir_count * 3
  iontorrent_evidence <- iontorrent_count * 2 + iontorrent_dir_count * 3
  
  # Find the platform with the most evidence
  platform_scores <- c(
    illumina = illumina_evidence,
    pacbio = pacbio_evidence,
    nanopore = nanopore_evidence,
    iontorrent = iontorrent_evidence
  )
  
  # Default to illumina if no strong evidence
  if (max(platform_scores) == 0) {
    return("illumina")
  }
  
  return(names(which.max(platform_scores)))
}

#' Get platform-specific patterns for file matching
#'
#' @param platform Character string indicating the sequencing platform
#' @param custom_patterns List of custom patterns (optional)
#'
#' @return List containing platform-specific patterns
#' @keywords internal
._getPlatformPatterns <- function(platform, custom_patterns = NULL) {
  # Default patterns
  default_patterns <- list(
    # Illumina patterns
    illumina = list(
      forward_patterns = c("_R1_", "_1\\.", "_READ1_", "_F\\.", "_forward\\."),
      reverse_patterns = c("_R2_", "_2\\.", "_READ2_", "_R\\.", "_reverse\\."),
      sample_patterns = c(
        "(.+?)(_S[0-9]+)?(_L[0-9]+)?(_R[12])",
        "(.+?)(_read[12])",
        "(.+?)(_[12])",
        "(.+?)_(forward|reverse)"
      )
    ),
    
    # PacBio patterns
    pacbio = list(
      forward_patterns = c("\\.subreads\\."),
      reverse_patterns = c("\\.scraps\\."),
      sample_patterns = c(
        "m([0-9]{5}_[0-9]{6})",
        "(.+?)(_subreads)"
      )
    ),
    
    # Nanopore patterns
    nanopore = list(
      forward_patterns = c("_1D_", "_pass_"),
      reverse_patterns = c("_2D_", "_fail_"),
      sample_patterns = c(
        "(.+?)(_pass)",
        "barcode([0-9]{2})",
        "(.+?)_(barcode[0-9]{2})"
      )
    ),
    
    # Ion Torrent patterns
    iontorrent = list(
      forward_patterns = c("_rawlib\\.", "_R1_", "_1\\."),
      reverse_patterns = c("_R2_", "_2\\."),
      sample_patterns = c(
        "IonXpress_([0-9]+)",
        "sample_([0-9]+)",
        "(.+?)(_rawlib)"
      )
    )
  )
  
  # Use custom patterns if provided, otherwise use defaults
  if (!is.null(custom_patterns)) {
    patterns <- custom_patterns
  } else {
    patterns <- default_patterns[[platform]]
    
    # Default to illumina patterns if platform not recognized
    if (is.null(patterns)) {
      warning("Platform '", platform, "' not recognized. Using Illumina patterns instead.")
      patterns <- default_patterns[["illumina"]]
    }
  }
  
  return(patterns)
}

#' Extract file information and metadata from file paths
#'
#' @param file_list Character vector of file paths
#' @param patterns List of platform-specific patterns
#' @param extract_metadata Logical indicating whether to extract metadata
#'
#' @return Tibble with file information and extracted metadata
#' @keywords internal
#' @importFrom dplyr mutate select %>%
#' @importFrom tibble tibble
#' @importFrom stringr str_extract str_match
#' @importFrom purrr map
._extractFileInfo <- function(file_list, patterns, extract_metadata) {
  # Create basic file info
  file_info <- tibble(
    file_path = file_list,
    file_name = basename(file_path),
    file_size = file.size(file_path),
    dir_path = dirname(file_path)
  )
  
  # Try each sample pattern until we get a match
  sample_names <- character(length(file_list))
  
  for (pattern in patterns$sample_patterns) {
    # For unmatched files, try the current pattern
    indices <- which(sample_names == "")
    if (length(indices) == 0) break
    
    matches <- str_match(file_info$file_name[indices], pattern)
    if (is.null(matches) || ncol(matches) < 2) next
    
    matched <- !is.na(matches[, 2])
    sample_names[indices[matched]] <- matches[matched, 2]
  }
  
  # Add sample names to the file info
  file_info$sample_name <- sample_names
  
  # Flag potential forward and reverse reads
  file_info <- file_info %>%
    mutate(
      potential_forward = ._detectPotentialFileType(file_name, patterns$forward_patterns),
      potential_reverse = ._detectPotentialFileType(file_name, patterns$reverse_patterns)
    )
  
  # Extract metadata if requested
  if (extract_metadata) {
    file_info <- file_info %>%
      mutate(metadata = map(file_name, ._extractMetadataFromFilename))
  } else {
    file_info <- file_info %>%
      mutate(metadata = map(seq_along(file_name), ~list()))
  }
  
  # Filter out files that couldn't be matched to a sample name
  file_info <- file_info %>%
    filter(sample_name != "")
  
  if (nrow(file_info) == 0) {
    warning("Could not extract sample names from any files using the platform-specific patterns")
  }
  
  return(file_info)
}

#' Detect whether a file matches forward or reverse patterns
#'
#' @param file_names Character vector of file names
#' @param patterns Character vector of patterns to match
#'
#' @return Logical vector indicating matches
#' @keywords internal
#' @importFrom stringr str_detect
._detectPotentialFileType <- function(file_names, patterns) {
  # Combine patterns with OR
  pattern <- paste(patterns, collapse = "|")
  str_detect(file_names, pattern)
}

#' Extract metadata from a filename using regular expressions
#'
#' @param filename Character string of filename
#'
#' @return List containing extracted metadata
#' @keywords internal
#' @importFrom stringr str_match str_extract_all
._extractMetadataFromFilename <- function(filename) {
  # Initialize empty metadata
  metadata <- list()
  
  # Extract potential date information (YYYYMMDD or YYYY-MM-DD formats)
  date_match <- str_match(filename, "([0-9]{4}[-_]?[0-9]{2}[-_]?[0-9]{2})")
  if (!is.na(date_match[1])) {
    metadata$date <- date_match[1]
  }
  
  # Extract potential lane information
  lane_match <- str_match(filename, "L([0-9]{3})")
  if (!is.na(lane_match[1])) {
    metadata$lane <- lane_match[2]
  }
  
  # Extract potential barcode information
  barcode_match <- str_match(filename, "(barcode[0-9]{2})")
  if (!is.na(barcode_match[1])) {
    metadata$barcode <- barcode_match[1]
  }
  
  # Extract sample number
  sample_match <- str_match(filename, "S([0-9]+)")
  if (!is.na(sample_match[1])) {
    metadata$sample_num <- sample_match[2]
  }
  
  # Extract flow cell ID (Illumina)
  flowcell_match <- str_match(filename, "([A-Z0-9]{9})")
  if (!is.na(flowcell_match[1])) {
    metadata$flowcell <- flowcell_match[1]
  }
  
  # Extract any numeric identifiers
  identifiers <- str_extract_all(filename, "([0-9]{3,})")[[1]]
  if (length(identifiers) > 0) {
    metadata$identifiers <- identifiers
  }
  
  return(metadata)
}

#' Pair FASTQ files with confidence scores using fuzzy matching
#'
#' @param file_info Tibble with file information
#' @param patterns List of platform-specific patterns
#' @param fuzzy_threshold Numeric threshold for fuzzy matching
#' @param platform Character string indicating the platform
#'
#' @return Tibble with paired files and confidence scores
#' @keywords internal
#' @importFrom dplyr filter select mutate left_join group_by ungroup %>%
#' @importFrom stringdist stringdist stringsim
#' @importFrom purrr map2_dbl
._pairFastqFiles <- function(file_info, patterns, fuzzy_threshold, platform) {
  # Separate potential forward and reverse files
  forward_files <- file_info %>%
    filter(potential_forward) %>%
    select(sample_name, forward = file_path, forward_size = file_size, 
           forward_metadata = metadata, forward_dir = dir_path)
  
  reverse_files <- file_info %>%
    filter(potential_reverse) %>%
    select(sample_name, reverse = file_path, reverse_size = file_size,
           reverse_metadata = metadata, reverse_dir = dir_path)
  
  # Group by sample name for initial matching
  potential_pairs <- forward_files %>%
    left_join(reverse_files, by = "sample_name")
  
  # Function to calculate confidence score for a potential pair
  calculate_confidence <- function(forward, reverse, forward_meta, reverse_meta, 
                                   forward_dir, reverse_dir, platform) {
    # Start with base confidence
    confidence <- 0.5
    
    # Get file base names
    forward_name <- basename(forward)
    reverse_name <- basename(reverse)
    
    # Calculate string similarity between filenames
    name_similarity <- stringdist::stringsim(forward_name, reverse_name, method = "jw")
    
    # Adjust for naming conventions based on platform
    if (platform == "illumina") {
      # Common Illumina pattern replacement (R1 -> R2)
      if (gsub("_R1_", "_R2_", forward_name) == reverse_name ||
          gsub("_1.", "_2.", forward_name) == reverse_name) {
        confidence <- confidence + 0.3
      }
    } else if (platform == "pacbio") {
      # PacBio specific adjustments
      if (gsub("\\.subreads\\.", ".scraps.", forward_name) == reverse_name) {
        confidence <- confidence + 0.3
      }
    } else if (platform == "nanopore") {
      # Nanopore specific adjustments
      if (gsub("_pass_", "_fail_", forward_name) == reverse_name) {
        confidence <- confidence + 0.3
      }
    }
    
    # Adjust based on name similarity
    confidence <- confidence + (name_similarity * 0.2)
    
    # Check if files are in the same directory
    if (forward_dir == reverse_dir) {
      confidence <- confidence + 0.1
    }
    
    # Check if file sizes are somewhat related (often reverse reads are smaller)
    if (!is.na(forward_size) && !is.na(reverse_size)) {
      size_ratio <- min(forward_size, reverse_size) / max(forward_size, reverse_size)
      if (size_ratio > 0.9) {
        confidence <- confidence + 0.1
      } else if (size_ratio > 0.7) {
        confidence <- confidence + 0.05
      }
    }
    
    # Check metadata consistency
    if (length(forward_meta) > 0 && length(reverse_meta) > 0) {
      common_fields <- intersect(names(forward_meta), names(reverse_meta))
      if (length(common_fields) > 0) {
        # Count matching fields
        matching_fields <- 0
        for (field in common_fields) {
          if (identical(forward_meta[[field]], reverse_meta[[field]])) {
            matching_fields <- matching_fields + 1
          }
        }
        
        # Adjust confidence based on metadata matches
        metadata_confidence <- matching_fields / max(1, length(common_fields)) * 0.1
        confidence <- confidence + metadata_confidence
      }
    }
    
    # Cap confidence at 1.0
    min(confidence, 1.0)
  }
  
  # Calculate confidence scores for all potential pairs
  paired_files <- potential_pairs %>%
    mutate(
      confidence = map2_dbl(
        forward, reverse, 
        ~calculate_confidence(.x, .y, forward_metadata, reverse_metadata, 
                             forward_dir, reverse_dir, platform)
      )
    ) %>%
    select(sample_name, forward, reverse, confidence, 
           forward_size, reverse_size, 
           metadata = forward_metadata, platform = forward_dir) %>%
    mutate(platform = platform)
  
  # Handle cases where no match was found based on sample name
  if (nrow(filter(paired_files, !is.na(reverse))) == 0) {
    # Try fuzzy matching across all forward and reverse files
    message("No exact matches found. Attempting fuzzy matching...")
    
    # Create all possible combinations for fuzzy matching
    cross_pairs <- expand.grid(
      forward_idx = seq_len(nrow(forward_files)),
      reverse_idx = seq_len(nrow(reverse_files)),
      stringsAsFactors = FALSE
    )
    
    # Calculate confidence for each potential cross-pair
    cross_confidence <- mapply(
      function(f_idx, r_idx) {
        calculate_confidence(
          forward_files$forward[f_idx],
          reverse_files$reverse[r_idx],
          forward_files$forward_metadata[f_idx],
          reverse_files$reverse_metadata[r_idx],
          forward_files$forward_dir[f_idx],
          reverse_files$reverse_dir[r_idx],
          platform
        )
      },
      cross_pairs$forward_idx,
      cross_pairs$reverse_idx
    )
    
    # Filter to keep only pairs above threshold
    valid_pairs <- which(cross_confidence >= fuzzy_threshold)
    
    if (length(valid_pairs) > 0) {
      # Create fuzzy-matched pairs
      fuzzy_pairs <- tibble(
        sample_name = paste0(
          forward_files$sample_name[cross_pairs$forward_idx[valid_pairs]],
          "_fuzzymatched"
        ),
        forward = forward_files$forward[cross_pairs$forward_idx[valid_pairs]],
        reverse = reverse_files$reverse[cross_pairs$reverse_idx[valid_pairs]],
        confidence = cross_confidence[valid_pairs],
        forward_size = forward_files$forward_size[cross_pairs$forward_idx[valid_pairs]],
        reverse_size = reverse_files$reverse_size[cross_pairs$reverse_idx[valid_pairs]],
        metadata = forward_files$forward_metadata[cross_pairs$forward_idx[valid_pairs]],
        platform = platform
      )
      
      # Combine with exact matches
      paired_files <- rbind(
        filter(paired_files, !is.na(reverse)),
        fuzzy_pairs
      )
    }
  }
  
  # Return paired files
  paired_files
}

#' Verify FASTQ content by sampling and analyzing reads
#'
#' @param file_data Tibble with file information
#' @param paired Logical indicating whether files are paired
#' @param verify_lines Integer number of lines to check
#'
#' @return The input tibble with added content verification columns
#' @keywords internal
#' @importFrom dplyr mutate %>%
#' @importFrom purrr map map_lgl
#' @importFrom readr read_lines
._verifyFastqContent <- function(file_data, paired, verify_lines) {
  # Function to verify FASTQ content and extract platform info
  verify_fastq <- function(file_path, verify_lines) {
    if (!._fileExists(file_path)) {
      return(list(valid = FALSE, avg_read_length = NA, header_format = NA))
    }
    
    tryCatch({
      # Read a sample of the file
      lines <- read_lines(file_path, n_max = verify_lines)
      
      if (length(lines) < 4) {
        return(list(valid = FALSE, avg_read_length = NA, header_format = NA))
      }
      
      # Check FASTQ format: header (@), sequence, separator (+), quality
      valid <- TRUE
      read_lengths <- numeric()
      headers <- character()
      
      for (i in seq(1, min(length(lines), verify_lines), by = 4)) {
        if (i + 3 > length(lines)) break
        
        if (!str_detect(lines[i], "^@") || !str_detect(lines[i+2], "^\\+")) {
          valid <- FALSE
          break
        }
        
        # Store read length
        read_lengths <- c(read_lengths, nchar(lines[i+1]))
        
        # Store header format
        headers <- c(headers, lines[i])
        
        # Check that sequence and quality lines have same length
        if (nchar(lines[i+1]) != nchar(lines[i+3])) {
          valid <- FALSE
          break
        }
      }
      
      # Extract header format
      header_format <- ._analyzeHeaderFormat(headers)
      
      return(list(
        valid = valid,
        avg_read_length = mean(read_lengths),
        header_format = header_format
      ))
    }, error = function(e) {
      message("Error validating file: ", file_path, " - ", e$message)
      return(list(valid = FALSE, avg_read_length = NA, header_format = NA))
    })
  }
  
  # Verify forward files
  if (paired) {
    verified_data <- file_data %>%
      mutate(
        forward_verification = map(forward, verify_fastq, verify_lines = verify_lines),
        forward_valid = map_lgl(forward_verification, "valid"),
        forward_avg_length = map_dbl(forward_verification, "avg_read_length"),
        forward_header = map_chr(forward_verification, "header_format")
      )
    
    # Verify reverse files
    verified_data <- verified_data %>%
      mutate(
        reverse_verification = map(reverse, verify_fastq, verify_lines = verify_lines),
        reverse_valid = map_lgl(reverse_verification, "valid"),
        reverse_avg_length = map_dbl(reverse_verification, "avg_read_length"),
        reverse_header = map_chr(reverse_verification, "header_format")
      )
    
    # Cross-validate header consistency
    verified_data <- verified_data %>%
      mutate(
        headers_consistent = forward_header == reverse_header,
        content_valid = forward_valid & reverse_valid
      )
  } else {
    # Single-end mode
    verified_data <- file_data %>%
      mutate(
        verification = map(file_path, verify_fastq, verify_lines = verify_lines),
        file_valid = map_lgl(verification, "valid"),
        avg_length = map_dbl(verification, "avg_read_length"),
        header_format = map_chr(verification, "header_format"),
        content_valid = file_valid
      )
  }
  
  # Return verified data
  return(verified_data)
}

#' Analyze FASTQ header format to determine platform characteristics
#'
#' @param headers Character vector of FASTQ headers
#'
#' @return Character string describing the header format
#' @keywords internal
#' @importFrom stringr str_detect
._analyzeHeaderFormat <- function(headers) {
  if (length(headers) == 0) {
    return("unknown")
  }
  
  # Sample a few headers for analysis
  sample_headers <- headers[1:min(10, length(headers))]
  
  # Illumina 1.8+ format: @EAS139:136:FC706VJ:2:2104:15343:197393
  if (any(str_detect(sample_headers, "^@[A-Z0-9]+:[0-9]+:[A-Z0-9]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+"))) {
    return("illumina_1.8+")
  }
  
  # Illumina 1.4+ format: @HWUSI-EAS100R:6:73:941:1973
  if (any(str_detect(sample_headers, "^@[A-Z0-9\\-]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+"))) {
    return("illumina_1.4+")
  }
  
  # PacBio format: @m54238_180625_052051/4194372/ccs
  if (any(str_detect(sample_headers, "^@m[0-9]+_[0-9]+_[0-9]+/[0-9]+/(ccs|subreads)"))) {
    return("pacbio")
  }
  
  # Nanopore format: @0c0aa36a-b0d2-4c89-9382-534507154868
  if (any(str_detect(sample_headers, "^@[a-z0-9]+\\-[a-z0-9]+\\-[a-z0-9]+\\-[a-z0-9]+\\-[a-z0-9]+"))) {
    return("nanopore")
  }
  
  # Ion Torrent format: @IARI0UB02G9O5A
  if (any(str_detect(sample_headers, "^@[A-Z0-9]{10,16}"))) {
    return("iontorrent")
  }
  
  return("unknown")
}

#' Cross-validate paired read assignments using Bayesian inference
#'
#' @param paired_files Tibble with paired file information
#' @param prior_weight Numeric weight to assign to priors
#'
#' @return The input tibble with updated confidence scores
#' @keywords internal
#' @importFrom dplyr mutate %>%
._bayesianValidatePairs <- function(paired_files, prior_weight = 0.3) {
  # Platform-specific priors (probability of correct pairing based on naming)
  platform_priors <- list(
    illumina = 0.95,
    pacbio = 0.9,
    nanopore = 0.85,
    iontorrent = 0.9,
    unknown = 0.7
  )
  
  # Get prior for this platform
  get_platform_prior <- function(platform) {
    if (platform %in% names(platform_priors)) {
      return(platform_priors[[platform]])
    }
    return(platform_priors[["unknown"]])
  }
  
  # Update confidence using Bayesian approach
  paired_files <- paired_files %>%
    mutate(
      platform_prior = map_dbl(platform, get_platform_prior),
      # Combine evidence with prior using weighted average
      bayesian_confidence = (confidence * (1 - prior_weight)) + (platform_prior * prior_weight)
    )
  
  # Replace original confidence with Bayesian confidence
  paired_files %>%
    mutate(confidence = bayesian_confidence) %>%
    select(-platform_prior, -bayesian_confidence)
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