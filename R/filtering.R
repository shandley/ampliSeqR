#' Filter and trim sequences based on quality metrics
#'
#' This function filters and trims sequences based on quality metrics using DADA2's
#' filterAndTrim function.
#'
#' @param fastq_files Tibble with FASTQ file information as produced by detectFastqFiles()
#' @param output_dir Character string specifying output directory for filtered files
#' @param truncLen Integer vector of length 2 giving truncation lengths for forward and reverse reads
#' @param maxEE Double vector of length 2 giving maximum expected errors for forward and reverse reads
#' @param truncQ Integer giving truncation quality threshold
#' @param maxN Integer giving maximum number of N bases allowed
#' @param rm.phix Logical indicating whether to remove PhiX reads
#' @param compress Logical indicating whether to gzip compress output files
#' @param multithread Logical or Integer indicating number of threads to use for processing
#' @param verbose Logical indicating whether to print verbose output
#'
#' @return A tibble with paths to filtered files and filtering statistics
#'
#' @export
#' @importFrom dada2 filterAndTrim
#' @importFrom dplyr %>% mutate bind_cols
#' @importFrom tibble tibble
#' @importFrom stringr str_replace
#'
#' @examples
#' \dontrun{
#' fastq_files <- detectFastqFiles("path/to/fastq/dir")
#' filtered_files <- filterAndTrimReads(fastq_files, "path/to/filtered/dir",
#'                                    truncLen = c(240, 200), maxEE = c(2, 2))
#' }
filterAndTrimReads <- function(fastq_files, 
                              output_dir,
                              truncLen = c(0, 0),
                              maxEE = c(2, 2),
                              truncQ = 2,
                              maxN = 0,
                              rm.phix = TRUE,
                              compress = TRUE,
                              multithread = FALSE,
                              verbose = TRUE) {
  
  # Input validation
  if (nrow(fastq_files) == 0) {
    stop("No FASTQ files provided")
  }
  
  if (!._dirExists(output_dir)) {
    if (verbose) message("Creating output directory: ", output_dir)
    if (!._createDirIfNotExists(output_dir)) {
      stop("Failed to create output directory: ", output_dir)
    }
  }
  
  # Check if we're working with paired-end data
  paired <- "reverse" %in% names(fastq_files) && all(!is.na(fastq_files$reverse))
  
  if (paired) {
    if (length(truncLen) < 2 || length(maxEE) < 2) {
      stop("For paired-end data, truncLen and maxEE must be vectors of length 2")
    }
  }
  
  # Create output file paths
  if (paired) {
    # Paired-end mode
    fastq_files <- fastq_files %>%
      mutate(
        filtered_forward = file.path(output_dir, paste0(sample_name, "_F_filtered.fastq.gz")),
        filtered_reverse = file.path(output_dir, paste0(sample_name, "_R_filtered.fastq.gz"))
      )
    
    # Set up filterAndTrim inputs  
    fwd_in <- fastq_files$forward
    rev_in <- fastq_files$reverse
    fwd_out <- fastq_files$filtered_forward
    rev_out <- fastq_files$filtered_reverse
    
    # Run filtering
    tryCatch({
      if (verbose) message("Filtering and trimming paired-end reads...")
      
      filter_stats <- filterAndTrim(
        fwd = fwd_in, fwd.out = fwd_out,
        rev = rev_in, rev.out = rev_out,
        truncLen = truncLen,
        maxEE = maxEE,
        truncQ = truncQ,
        maxN = maxN,
        rm.phix = rm.phix,
        compress = compress,
        multithread = multithread,
        verbose = verbose
      )
      
      # Convert filter_stats to tibble and add to results
      filter_stats_tbl <- as_tibble(filter_stats, rownames = "file") %>%
        mutate(file = basename(file)) %>%
        rename(reads_in = reads.in, reads_out = reads.out)
      
      # Split into forward and reverse stats
      fwd_stats <- filter_stats_tbl %>%
        filter(str_detect(file, "_F_filtered")) %>%
        rename(forward_reads_in = reads_in, forward_reads_out = reads_out) %>%
        mutate(sample_name = str_replace(file, "_F_filtered.fastq.gz", ""))
      
      rev_stats <- filter_stats_tbl %>%
        filter(str_detect(file, "_R_filtered")) %>%
        rename(reverse_reads_in = reads_in, reverse_reads_out = reads_out) %>%
        mutate(sample_name = str_replace(file, "_R_filtered.fastq.gz", ""))
      
      # Combine stats
      stats_combined <- left_join(fwd_stats, rev_stats, by = "sample_name") %>%
        select(-file.x, -file.y)
      
      # Join stats with files
      result <- left_join(fastq_files, stats_combined, by = "sample_name")
      
    }, error = function(e) {
      stop("Error during filtering and trimming: ", e$message)
    })
    
  } else {
    # Single-end mode
    fastq_files <- fastq_files %>%
      mutate(
        filtered_forward = file.path(output_dir, paste0(sample_name, "_filtered.fastq.gz"))
      )
    
    # Set up filterAndTrim inputs  
    fwd_in <- fastq_files$forward
    fwd_out <- fastq_files$filtered_forward
    
    # Run filtering
    tryCatch({
      if (verbose) message("Filtering and trimming single-end reads...")
      
      filter_stats <- filterAndTrim(
        fwd = fwd_in, fwd.out = fwd_out,
        truncLen = truncLen[1],
        maxEE = maxEE[1],
        truncQ = truncQ,
        maxN = maxN,
        rm.phix = rm.phix,
        compress = compress,
        multithread = multithread,
        verbose = verbose
      )
      
      # Convert filter_stats to tibble and add to results
      filter_stats_tbl <- as_tibble(filter_stats, rownames = "file") %>%
        mutate(file = basename(file)) %>%
        rename(forward_reads_in = reads.in, forward_reads_out = reads.out) %>%
        mutate(sample_name = str_replace(file, "_filtered.fastq.gz", ""))
      
      # Join stats with files
      result <- left_join(fastq_files, filter_stats_tbl, by = "sample_name") %>%
        select(-file)
      
    }, error = function(e) {
      stop("Error during filtering and trimming: ", e$message)
    })
  }
  
  if (verbose) {
    total_in <- sum(result$forward_reads_in, na.rm = TRUE)
    total_out <- sum(result$forward_reads_out, na.rm = TRUE)
    retention <- round(100 * total_out / total_in, 1)
    
    message("Filtered ", nrow(result), " samples (", 
            total_in, " reads in, ", total_out, " reads out, ",
            retention, "% retention)")
  }
  
  return(result)
}

#' Optimize filtering parameters
#'
#' This function helps identify optimal filtering parameters by testing different
#' combinations of truncation lengths and maxEE values.
#'
#' @param fastq_files Tibble with FASTQ file information as produced by detectFastqFiles()
#' @param output_dir Character string specifying temporary output directory for filtered files
#' @param truncLen_range List of integer vectors to test for truncation lengths
#' @param maxEE_range List of double vectors to test for maximum expected errors
#' @param sample_size Integer indicating number of samples to use for optimization
#' @param paired Logical indicating whether files are paired
#'
#' @return A tibble with filtering results for different parameter combinations
#'
#' @export
#' @importFrom dplyr %>% group_by summarize
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' fastq_files <- detectFastqFiles("path/to/fastq/dir")
#' truncLen_range <- list(c(240, 200), c(220, 180), c(200, 160))
#' maxEE_range <- list(c(2, 2), c(2, 5), c(5, 5))
#' opt_params <- optimizeFilteringParams(fastq_files, "path/to/temp/dir",
#'                                      truncLen_range, maxEE_range)
#' }
optimizeFilteringParams <- function(fastq_files, 
                                   output_dir,
                                   truncLen_range = list(c(240, 200), c(220, 180), c(200, 160)),
                                   maxEE_range = list(c(2, 2), c(2, 5), c(5, 5)),
                                   sample_size = 2,
                                   paired = TRUE) {
  
  # Input validation
  if (nrow(fastq_files) == 0) {
    stop("No FASTQ files provided")
  }
  
  if (sample_size > nrow(fastq_files)) {
    sample_size <- nrow(fastq_files)
    warning("sample_size greater than number of samples, using all samples")
  }
  
  # Create temporary directory
  temp_dir <- file.path(output_dir, "param_opt")
  if (!._createDirIfNotExists(temp_dir)) {
    stop("Failed to create temporary directory: ", temp_dir)
  }
  
  # Sample files for testing
  set.seed(42)  # For reproducibility
  sample_idx <- sample(1:nrow(fastq_files), sample_size)
  sample_files <- fastq_files[sample_idx, ]
  
  # Function to test one parameter combination
  ._testParams <- function(truncLen, maxEE) {
    param_dir <- file.path(temp_dir, 
                          paste0("truncLen_", paste(truncLen, collapse = "_"), 
                                 "_maxEE_", paste(maxEE, collapse = "_")))
    
    if (!._createDirIfNotExists(param_dir)) {
      stop("Failed to create parameter directory: ", param_dir)
    }
    
    message("Testing parameters: truncLen = ", paste(truncLen, collapse = ", "), 
            ", maxEE = ", paste(maxEE, collapse = ", "))
    
    # Run filtering with these parameters
    filter_results <- tryCatch({
      filterAndTrimReads(
        sample_files, 
        param_dir,
        truncLen = truncLen,
        maxEE = maxEE,
        verbose = FALSE
      )
    }, error = function(e) {
      warning("Error with parameters (truncLen = ", paste(truncLen, collapse = ", "), 
              ", maxEE = ", paste(maxEE, collapse = ", "), "): ", e$message)
      return(NULL)
    })
    
    if (is.null(filter_results)) {
      return(tibble(
        truncLen_fwd = truncLen[1],
        truncLen_rev = if(paired) truncLen[2] else NA_integer_,
        maxEE_fwd = maxEE[1],
        maxEE_rev = if(paired) maxEE[2] else NA_real_,
        total_reads_in = NA_integer_,
        total_reads_out = NA_integer_,
        retention_pct = NA_real_,
        error = TRUE
      ))
    } else {
      # Calculate retention statistics
      total_in <- sum(filter_results$forward_reads_in, na.rm = TRUE)
      total_out <- sum(filter_results$forward_reads_out, na.rm = TRUE)
      retention <- if(total_in > 0) 100 * total_out / total_in else 0
      
      return(tibble(
        truncLen_fwd = truncLen[1],
        truncLen_rev = if(paired) truncLen[2] else NA_integer_,
        maxEE_fwd = maxEE[1],
        maxEE_rev = if(paired) maxEE[2] else NA_real_,
        total_reads_in = total_in,
        total_reads_out = total_out,
        retention_pct = retention,
        error = FALSE
      ))
    }
  }
  
  # Test all parameter combinations
  results <- tibble(
    truncLen = truncLen_range,
    maxEE = rep(list(NA), length(truncLen_range))
  )
  
  all_results <- vector("list", length(truncLen_range) * length(maxEE_range))
  counter <- 1
  
  for (i in seq_along(truncLen_range)) {
    for (j in seq_along(maxEE_range)) {
      all_results[[counter]] <- ._testParams(truncLen_range[[i]], maxEE_range[[j]])
      counter <- counter + 1
    }
  }
  
  combined_results <- bind_rows(all_results)
  
  # Clean up temporary files if requested
  # Uncomment to enable cleanup
  # unlink(temp_dir, recursive = TRUE)
  
  # Sort by retention and return
  optimal_params <- combined_results %>%
    filter(!error) %>%
    arrange(desc(retention_pct))
  
  return(optimal_params)
}