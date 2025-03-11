#!/usr/bin/env Rscript

# Example script for customizing sequence length parameters
# This addresses the "no sequences retained" error

library(ampliSeqR)

# Example code to analyze your data with custom sequence length parameters
analyze_data <- function() {
  # Your existing code to detect and filter FASTQ files
  # fastq_files <- detectFastqFiles("path/to/fastq/dir")
  # filtered_files <- filterAndTrimReads(fastq_files, "path/to/filtered/dir")
  
  # Determine your expected amplicon length range
  # This depends on your specific primers and 16S region!
  #
  # For V3-V4 region: ~450bp
  # For V4 region: ~250bp
  # For V1-V3 region: ~500bp
  # For full-length 16S: ~1500bp
  
  # Example for V4 region (default in package)
  min_length <- 240
  max_length <- 260
  
  # Example for V3-V4 region
  # min_length <- 430
  # max_length <- 470
  
  # Run the DADA2 pipeline with your custom length parameters
  # results <- runDADA2Pipeline(
  #   filtered_files,
  #   n_reads = 1e8,
  #   min_length = min_length,
  #   max_length = max_length
  # )
  
  cat("To fix the 'No sequences remained after length filtering' error:\n\n")
  cat("1. Identify the expected amplicon length for your primers/protocol\n")
  cat("2. Set min_length and max_length to span this range with some buffer\n")
  cat("3. For example, for V3-V4, you might use min_length=430, max_length=470\n\n")
  cat("Current error shows filtering with range 380-430 bp, which isn't matching your amplicons.\n")
  cat("Try customizing these parameters in your runDADA2Pipeline call.\n")
}

analyze_data()