#' Learn error rates using DADA2
#'
#' This function learns error rates from filtered sequence files using DADA2's
#' learnErrors function.
#'
#' @param filtered_files Tibble with filtered file information as produced by filterAndTrimReads()
#' @param n_reads Number of bases to use for error learning (default: 1e8)
#' @param n_iters Number of iterations for error learning (default: 12)
#' @param multithread Logical or Integer indicating number of threads to use (default: FALSE)
#' @param verbose Logical indicating whether to print verbose output (default: TRUE)
#' @param paired Logical indicating whether data is paired-end (default: TRUE)
#'
#' @return A list with error models for forward and reverse reads
#'
#' @export
#' @importFrom dada2 learnErrors
#' @importFrom dplyr filter
#'
#' @examples
#' \dontrun{
#' fastq_files <- detectFastqFiles("path/to/fastq/dir")
#' filtered_files <- filterAndTrimReads(fastq_files, "path/to/filtered/dir")
#' error_models <- learnErrors(filtered_files)
#' }
learnErrorRates <- function(filtered_files,
                           n_reads = 1e8,
                           n_iters = 12,
                           multithread = FALSE,
                           verbose = TRUE,
                           paired = TRUE) {
  
  # Input validation
  if (nrow(filtered_files) == 0) {
    stop("No filtered files provided")
  }
  
  # Check if we're working with paired-end data
  paired <- paired && "filtered_reverse" %in% names(filtered_files) && 
    all(!is.na(filtered_files$filtered_reverse))
  
  # Learn errors for forward reads
  if (verbose) message("Learning error rates for forward reads...")
  err_f <- tryCatch({
    learnErrors(filtered_files$filtered_forward, 
               multithread = multithread,
               nbases = n_reads,
               MAX_CONSIST = n_iters,
               verbose = verbose)
  }, error = function(e) {
    stop("Error learning forward read error rates: ", e$message)
  })
  
  # Learn errors for reverse reads if paired
  if (paired) {
    if (verbose) message("Learning error rates for reverse reads...")
    err_r <- tryCatch({
      learnErrors(filtered_files$filtered_reverse, 
                 multithread = multithread,
                 nbases = n_reads,
                 MAX_CONSIST = n_iters,
                 verbose = verbose)
    }, error = function(e) {
      stop("Error learning reverse read error rates: ", e$message)
    })
    
    return(list(forward = err_f, reverse = err_r))
  } else {
    return(list(forward = err_f))
  }
}

#' Infer sample composition using DADA2
#'
#' This function infers the sample composition from filtered sequence files using DADA2's
#' dada function and the error models learned with learnErrorRates.
#'
#' @param filtered_files Tibble with filtered file information as produced by filterAndTrimReads()
#' @param error_models List with error models as produced by learnErrorRates()
#' @param pool Logical or character indicating whether to pool sequences (default: FALSE)
#' @param multithread Logical or Integer indicating number of threads to use (default: FALSE)
#' @param verbose Logical indicating whether to print verbose output (default: TRUE)
#' @param paired Logical indicating whether data is paired-end (default: TRUE)
#'
#' @return A list with dada objects for each sample
#'
#' @export
#' @importFrom dada2 dada mergePairs
#' @importFrom dplyr %>% mutate
#' @importFrom purrr map
#'
#' @examples
#' \dontrun{
#' fastq_files <- detectFastqFiles("path/to/fastq/dir")
#' filtered_files <- filterAndTrimReads(fastq_files, "path/to/filtered/dir")
#' error_models <- learnErrorRates(filtered_files)
#' dada_objects <- inferSampleComposition(filtered_files, error_models)
#' }
inferSampleComposition <- function(filtered_files,
                                  error_models,
                                  pool = FALSE,
                                  multithread = FALSE,
                                  verbose = TRUE,
                                  paired = TRUE) {
  
  # Input validation
  if (nrow(filtered_files) == 0) {
    stop("No filtered files provided")
  }
  
  if (is.null(error_models$forward)) {
    stop("Forward error model not found in error_models")
  }
  
  # Check if we're working with paired-end data
  paired <- paired && "filtered_reverse" %in% names(filtered_files) && 
    all(!is.na(filtered_files$filtered_reverse)) &&
    !is.null(error_models$reverse)
  
  if (paired && is.null(error_models$reverse)) {
    stop("Paired-end mode requested but reverse error model not found")
  }
  
  # Process forward reads
  if (verbose) message("Denoising forward reads...")
  dada_fwd <- dada(filtered_files$filtered_forward, 
                  err = error_models$forward, 
                  pool = pool,
                  multithread = multithread)
  
  # Process reverse reads and merge if paired
  if (paired) {
    if (verbose) message("Denoising reverse reads...")
    dada_rev <- dada(filtered_files$filtered_reverse, 
                    err = error_models$reverse, 
                    pool = pool,
                    multithread = multithread)
    
    if (verbose) message("Merging paired reads...")
    mergers <- mergePairs(dada_fwd, filtered_files$filtered_forward,
                         dada_rev, filtered_files$filtered_reverse,
                         verbose = verbose)
    
    return(list(
      forward = dada_fwd,
      reverse = dada_rev,
      merged = mergers
    ))
  } else {
    return(list(
      forward = dada_fwd
    ))
  }
}

#' Generate ASV table from DADA2 results
#'
#' This function generates an ASV (Amplicon Sequence Variant) table from the
#' results of DADA2 processing.
#'
#' @param dada_results List with dada objects as produced by inferSampleComposition()
#' @param filtered_files Tibble with filtered file information
#' @param min_length Minimum acceptable sequence length (default: 0)
#' @param max_length Maximum acceptable sequence length (default: Inf)
#' @param remove_chimeras Logical indicating whether to remove chimeras (default: TRUE)
#' @param method Character string specifying chimera removal method (default: "consensus")
#' @param verbose Logical indicating whether to print verbose output (default: TRUE)
#' @param paired Logical indicating whether data is paired-end (default: TRUE)
#'
#' @return A list containing the ASV table, ASV sequences, and summary statistics
#'
#' @export
#' @importFrom dada2 makeSequenceTable removeBimeraDenovo
#' @importFrom Biostrings DNAStringSet writeXStringSet
#' @importFrom dplyr %>% mutate
#'
#' @examples
#' \dontrun{
#' fastq_files <- detectFastqFiles("path/to/fastq/dir")
#' filtered_files <- filterAndTrimReads(fastq_files, "path/to/filtered/dir")
#' error_models <- learnErrorRates(filtered_files)
#' dada_objects <- inferSampleComposition(filtered_files, error_models)
#' asv_results <- generateASVTable(dada_objects, filtered_files)
#' }
generateASVTable <- function(dada_results,
                            filtered_files,
                            min_length = 0,
                            max_length = Inf,
                            remove_chimeras = TRUE,
                            method = "consensus",
                            verbose = TRUE,
                            paired = TRUE) {
  
  # Input validation
  if (is.null(dada_results$forward)) {
    stop("Forward dada results not found")
  }
  
  # Check data type (paired or single)
  paired <- paired && !is.null(dada_results$merged)
  
  # Make sequence table
  if (paired) {
    if (verbose) message("Creating sequence table from merged reads...")
    seqtab <- makeSequenceTable(dada_results$merged)
  } else {
    if (verbose) message("Creating sequence table from forward reads...")
    seqtab <- makeSequenceTable(dada_results$forward)
  }
  
  # Filter by sequence length
  if (min_length > 0 || max_length < Inf) {
    seq_lengths <- nchar(colnames(seqtab))
    filtered_seqs <- seq_lengths >= min_length & seq_lengths <= max_length
    n_retained <- sum(filtered_seqs)
    
    if (verbose) {
      message("Filtered sequences by length (", min_length, "-", max_length, " bp): ",
              n_retained, " of ", length(seq_lengths), " sequences retained")
    }
    
    # Stop if no sequences were retained
    if (n_retained == 0) {
      stop("No sequences remained after length filtering. Try adjusting min_length and max_length parameters to match your amplicon size.")
    }
    
    seqtab <- seqtab[, filtered_seqs]
  }
  
  # Remove chimeras
  if (remove_chimeras) {
    if (verbose) message("Removing chimeric sequences...")
    
    original_cols <- ncol(seqtab)
    seqtab.nochim <- removeBimeraDenovo(seqtab, method = method, multithread = TRUE, verbose = verbose)
    
    if (verbose) {
      chimera_cols <- original_cols - ncol(seqtab.nochim)
      message("Removed ", chimera_cols, " chimeric sequences (", 
              round(chimera_cols/original_cols*100, 1), "%)")
    }
    
    seqtab <- seqtab.nochim
  }
  
  # Create DNA string set
  asv_sequences <- DNAStringSet(colnames(seqtab))
  names(asv_sequences) <- paste0("ASV", seq_along(asv_sequences))
  
  # Rename ASVs in sequence table
  seqtab_renamed <- seqtab
  colnames(seqtab_renamed) <- names(asv_sequences)
  
  # Calculate summary statistics
  sample_reads <- rowSums(seqtab_renamed)
  total_reads <- sum(sample_reads)
  total_asvs <- ncol(seqtab_renamed)
  
  # Create ASV prevalence and abundance stats
  asv_prevalence <- colSums(seqtab_renamed > 0)
  asv_abundance <- colSums(seqtab_renamed)
  
  asv_stats <- tibble(
    asv_id = names(asv_sequences),
    sequence = as.character(asv_sequences),
    length = nchar(sequence),
    abundance = asv_abundance,
    prevalence = asv_prevalence,
    rel_abundance = abundance / total_reads
  )
  
  # Extract simplified sample names (remove file extension if present)
  sample_names <- rownames(seqtab_renamed)
  # If the sample names include file extensions, clean them up
  sample_names_clean <- gsub("_F_filtered\\.fastq\\.gz$", "", sample_names)
  sample_names_clean <- gsub("_filtered\\.fastq\\.gz$", "", sample_names_clean)
  sample_names_clean <- gsub("\\.fastq\\.gz$", "", sample_names_clean)
  sample_names_clean <- gsub("\\.fastq$", "", sample_names_clean)
  sample_names_clean <- basename(sample_names_clean)  # In case there are paths
  
  # Return results
  if (verbose) {
    message("Final ASV table: ", nrow(seqtab_renamed), " samples x ", ncol(seqtab_renamed), " ASVs")
    message("Total reads in ASV table: ", total_reads)
  }
  
  return(list(
    asv_table = seqtab_renamed,
    asv_sequences = asv_sequences,
    asv_stats = asv_stats,
    sample_stats = tibble(
      sample = sample_names_clean,  # Use cleaned sample names
      original_sample = sample_names,  # Keep original for reference
      total_reads = sample_reads
    )
  ))
}

#' Run complete DADA2 pipeline
#'
#' This function runs the complete DADA2 pipeline from filtered files to ASV table.
#'
#' @param filtered_files Tibble with filtered file information as produced by filterAndTrimReads()
#' @param n_reads Number of bases to use for error learning (default: 1e8)
#' @param pool Logical or character indicating whether to pool sequences (default: FALSE)
#' @param min_length Minimum acceptable sequence length (default: 240)
#' @param max_length Maximum acceptable sequence length (default: 260)
#' @param remove_chimeras Logical indicating whether to remove chimeras (default: TRUE)
#' @param multithread Logical or Integer indicating number of threads to use (default: FALSE)
#' @param verbose Logical indicating whether to print verbose output (default: TRUE)
#' @param paired Logical indicating whether data is paired-end (default: TRUE)
#'
#' @return A list containing the ASV table, ASV sequences, and summary statistics
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fastq_files <- detectFastqFiles("path/to/fastq/dir")
#' filtered_files <- filterAndTrimReads(fastq_files, "path/to/filtered/dir")
#' results <- runDADA2Pipeline(filtered_files)
#' }
runDADA2Pipeline <- function(filtered_files,
                            n_reads = 1e8,
                            pool = FALSE,
                            min_length = 240,
                            max_length = 260,
                            remove_chimeras = TRUE,
                            multithread = FALSE,
                            verbose = TRUE,
                            paired = TRUE) {
  
  # Input validation
  if (nrow(filtered_files) == 0) {
    stop("No filtered files provided")
  }
  
  # Check for paired-end files
  paired <- paired && "filtered_reverse" %in% names(filtered_files) && 
    all(!is.na(filtered_files$filtered_reverse))
  
  # Step 1: Learn error rates
  if (verbose) message("Step 1: Learning error rates...")
  error_models <- learnErrorRates(
    filtered_files = filtered_files,
    n_reads = n_reads,
    multithread = multithread,
    verbose = verbose,
    paired = paired
  )
  
  # Step 2: Infer sample composition
  if (verbose) message("Step 2: Inferring sample composition...")
  dada_results <- inferSampleComposition(
    filtered_files = filtered_files,
    error_models = error_models,
    pool = pool,
    multithread = multithread,
    verbose = verbose,
    paired = paired
  )
  
  # Step 3: Generate ASV table
  if (verbose) message("Step 3: Generating ASV table...")
  asv_results <- generateASVTable(
    dada_results = dada_results,
    filtered_files = filtered_files,
    min_length = min_length,
    max_length = max_length,
    remove_chimeras = remove_chimeras,
    verbose = verbose,
    paired = paired
  )
  
  # Add error models and dada objects to results
  asv_results$error_models <- error_models
  asv_results$dada_results <- dada_results
  
  return(asv_results)
}