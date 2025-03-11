#' Calculate quality statistics for FASTQ files
#'
#' This function analyzes FASTQ file quality and returns quality metrics
#' for each sample, including per-base quality scores and read length distribution.
#'
#' @param fastq_files Tibble with FASTQ file information as produced by detectFastqFiles()
#' @param n_reads Number of reads to sample for quality analysis (default: 10000)
#' @param paired Logical indicating whether files are paired (default: TRUE)
#'
#' @return A list containing quality objects for forward and reverse reads (if paired)
#'
#' @export
#' @importFrom ShortRead qa FastqQuality
#' @importFrom dplyr %>% mutate
#' @importFrom purrr map
#'
#' @examples
#' \dontrun{
#' fastq_files <- detectFastqFiles("path/to/fastq/dir")
#' quality_data <- calculateQualityStats(fastq_files)
#' }
calculateQualityStats <- function(fastq_files, n_reads = 10000, paired = TRUE) {
  
  # Input validation
  if (nrow(fastq_files) == 0) {
    stop("No FASTQ files provided")
  }
  
  if (paired && !"reverse" %in% names(fastq_files)) {
    warning("Paired mode requested but reverse files not found in input")
    paired <- FALSE
  }
  
  # Function to calculate quality metrics for a file
  ._calculateFileQuality <- function(file_path, n_reads) {
    message("Analyzing quality for: ", basename(file_path))
    
    tryCatch({
      # Use ShortRead for quality analysis
      qa_data <- qa(file_path, n = n_reads)
      return(qa_data)
    }, error = function(e) {
      warning("Failed to analyze quality for file: ", basename(file_path), " - ", e$message)
      return(NULL)
    })
  }
  
  # Calculate quality for forward reads
  forward_quality <- fastq_files %>%
    filter(!is.na(forward)) %>%
    mutate(quality = map(forward, ._calculateFileQuality, n_reads = n_reads))
  
  # Calculate quality for reverse reads if paired
  if (paired) {
    reverse_quality <- fastq_files %>%
      filter(!is.na(reverse)) %>%
      mutate(quality = map(reverse, ._calculateFileQuality, n_reads = n_reads))
    
    quality_data <- list(
      forward = forward_quality,
      reverse = reverse_quality
    )
  } else {
    quality_data <- list(
      forward = forward_quality
    )
  }
  
  return(quality_data)
}

#' Analyze read length distribution
#'
#' This function analyzes the distribution of read lengths in FASTQ files.
#'
#' @param fastq_files Tibble with FASTQ file information as produced by detectFastqFiles()
#' @param n_reads Number of reads to sample for length analysis (default: 10000)
#' @param paired Logical indicating whether files are paired (default: TRUE)
#'
#' @return A tibble with read length distributions for each sample
#'
#' @export
#' @importFrom ShortRead readFastq sread
#' @importFrom Biostrings width
#' @importFrom tibble tibble
#' @importFrom dplyr %>% mutate group_by summarize
#' @importFrom purrr map
#'
#' @examples
#' \dontrun{
#' fastq_files <- detectFastqFiles("path/to/fastq/dir")
#' length_data <- analyzeReadLengths(fastq_files)
#' }
analyzeReadLengths <- function(fastq_files, n_reads = 10000, paired = TRUE) {
  
  # Input validation
  if (nrow(fastq_files) == 0) {
    stop("No FASTQ files provided")
  }
  
  # Function to analyze read lengths
  ._getReadLengths <- function(file_path, n_reads) {
    message("Analyzing lengths for: ", basename(file_path))
    
    tryCatch({
      # Sample reads
      reads <- readFastq(file_path, n = n_reads)
      
      # Get read lengths
      lengths <- width(sread(reads))
      
      # Create distribution
      length_dist <- tibble(
        length = lengths
      ) %>%
        group_by(length) %>%
        summarize(count = n(), .groups = "drop") %>%
        mutate(proportion = count / sum(count))
      
      return(length_dist)
    }, error = function(e) {
      warning("Failed to analyze lengths for file: ", basename(file_path), " - ", e$message)
      return(tibble(length = numeric(0), count = numeric(0), proportion = numeric(0)))
    })
  }
  
  # Analyze forward reads
  length_data <- fastq_files %>%
    filter(!is.na(forward)) %>%
    mutate(
      sample_name = sample_name,
      read_type = "forward",
      length_dist = map(forward, ._getReadLengths, n_reads = n_reads)
    )
  
  # Analyze reverse reads if paired
  if (paired && "reverse" %in% names(fastq_files)) {
    reverse_length_data <- fastq_files %>%
      filter(!is.na(reverse)) %>%
      mutate(
        sample_name = sample_name,
        read_type = "reverse",
        length_dist = map(reverse, ._getReadLengths, n_reads = n_reads)
      )
    
    length_data <- bind_rows(length_data, reverse_length_data)
  }
  
  return(length_data)
}

#' Calculate summary statistics from quality data
#'
#' @param quality_data Quality data as returned by calculateQualityStats()
#'
#' @return Tibble with summary quality statistics for each sample and read type
#'
#' @export
#' @importFrom dplyr mutate summarize group_by
#' @importFrom purrr map map_dbl
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' fastq_files <- detectFastqFiles("path/to/fastq/dir")
#' quality_data <- calculateQualityStats(fastq_files)
#' quality_summary <- summarizeQualityStats(quality_data)
#' }
summarizeQualityStats <- function(quality_data) {
  
  # Extract summary stats from ShortRead qa objects
  ._extractQualitySummary <- function(qa_obj) {
    if (is.null(qa_obj)) {
      return(tibble(
        mean_quality = NA_real_,
        median_quality = NA_real_,
        min_quality = NA_real_,
        max_quality = NA_real_
      ))
    }
    
    # Get per-cycle quality scores
    per_cycle <- qa_obj[["perCycle"]]$quality
    
    # Calculate summary stats for each position
    pos_summary <- per_cycle %>%
      group_by(Cycle) %>%
      summarize(
        mean_quality = mean(Score),
        median_quality = median(Score),
        min_quality = min(Score),
        max_quality = max(Score),
        .groups = "drop"
      )
    
    # Calculate overall stats
    overall_summary <- tibble(
      mean_quality = mean(per_cycle$Score),
      median_quality = median(per_cycle$Score),
      min_quality = min(per_cycle$Score),
      max_quality = max(per_cycle$Score)
    )
    
    return(list(
      overall = overall_summary,
      per_position = pos_summary
    ))
  }
  
  # Process forward reads
  forward_summary <- quality_data$forward %>%
    mutate(
      read_type = "forward",
      quality_summary = map(quality, ._extractQualitySummary)
    )
  
  # Process reverse reads if available
  if ("reverse" %in% names(quality_data)) {
    reverse_summary <- quality_data$reverse %>%
      mutate(
        read_type = "reverse",
        quality_summary = map(quality, ._extractQualitySummary)
      )
    
    quality_summary <- bind_rows(forward_summary, reverse_summary)
  } else {
    quality_summary <- forward_summary
  }
  
  return(quality_summary)
}