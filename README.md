## 🚧 Under Construction 🚧

This project is currently under active development. Some features may be incomplete or subject to change.

# ampliSeqR

An R package for efficient 16S amplicon sequencing analysis with hardware acceleration.

## Overview

ampliSeqR provides tools for analyzing 16S rRNA gene amplicon sequencing data, with a focus on performance and hardware acceleration. This package wraps and extends the functionality of DADA2 to simplify and optimize the analysis workflow.

## Features

- Advanced file detection and organization with fuzzy matching capabilities
- Support for multiple sequencing platforms (Illumina, PacBio, Nanopore, Ion Torrent)
- Bayesian confidence scoring for file pairing
- Directory structure analysis for multi-batch experiments
- Quality control and visualization of sequencing data
- Optimized filtering and trimming of sequences
- Accelerated DADA2 implementation for error learning and sample inference
- ASV table generation and statistics
- Visualization functions for quality assessment and results

## Installation

```r
# Install dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("dada2", "Biostrings", "ShortRead"))

# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# Install the package
# Replace with your actual repository details when available
devtools::install_github("username/ampliSeqR")
```

## Quick Start

```r
library(ampliSeqR)
library(ggplot2)  # For visualizations

# Step 1: Detect and validate FASTQ files
fastq_dir <- "path/to/fastq/files"  # Replace with your actual path

# Choose between standard or advanced file detection:
# Standard detection for simple naming conventions:
fastq_files <- detectFastqFiles(fastq_dir, paired = TRUE)

# OR use advanced detection for complex file naming or multi-platform studies:
fastq_files <- detectFastqFilesAdvanced(
  fastq_dir, 
  paired = TRUE,
  platform = "auto",      # Auto-detect the sequencing platform
  recursive = TRUE,       # Search in nested directories
  fuzzy_threshold = 0.2,  # Allow fuzzy matching for filenames
  extract_metadata = TRUE # Extract metadata from filenames
)

validated_files <- validateFastqFiles(fastq_files)

# Step 2: Analyze quality
quality_data <- calculateQualityStats(validated_files)
# View quality profiles
quality_plot <- plotQualityProfiles(quality_data)
print(quality_plot)

# Analyze read lengths
length_data <- analyzeReadLengths(validated_files)
length_plot <- plotReadLengths(length_data)
print(length_plot)

# Step 3: Filter and trim reads
# Create output directory for filtered files
filtered_dir <- "path/to/filtered/files"  # Replace with your actual path
dir.create(filtered_dir, showWarnings = FALSE, recursive = TRUE)

# Optional: Optimize filtering parameters 
# (you can skip this and use default parameters)
truncLen_range <- list(c(240, 200), c(220, 180), c(200, 160))
maxEE_range <- list(c(2, 2), c(2, 5), c(5, 5))
opt_params <- optimizeFilteringParams(
  validated_files, 
  filtered_dir, 
  truncLen_range = truncLen_range,
  maxEE_range = maxEE_range,
  sample_size = 2  # Use 2 samples for optimization
)

# Use optimal parameters or default ones
filtered_files <- filterAndTrimReads(
  validated_files, 
  filtered_dir,
  truncLen = c(240, 200),  # Adjust based on your data
  maxEE = c(2, 2)
)

# Step 4: Run DADA2 pipeline (this performs the core analysis)
results <- runDADA2Pipeline(
  filtered_files,
  nbases = 1e8,  # Number of bases to process
  pool = FALSE,
  min_length = 240,  # Adjust based on your expected amplicon size
  max_length = 260,
  remove_chimeras = TRUE,
  multithread = TRUE
)

# Step 5: Visualize results
# Plot ASV abundances
asv_plot <- plotASVAbundances(results, top_n = 15, log_scale = TRUE)
print(asv_plot)

# Skip directly to saving results for downstream analysis

# Step 6: Save results for downstream analysis
# Save ASV table
write.csv(results$asv_table, "asv_table.csv")

# Save ASV sequences
Biostrings::writeXStringSet(
  results$asv_sequences, 
  "asv_sequences.fasta"
)

# Save ASV statistics
write.csv(results$asv_stats, "asv_stats.csv")
```

## Complete Workflow Example

For a more detailed example, see the included demo script:

```r
# From R console
file.edit(system.file("demo_script.R", package = "ampliSeqR"))
```

## Hardware Acceleration

ampliSeqR provides automatic hardware detection and optimization to accelerate analyses across different computing environments:

```r
# Detect available hardware resources
hw_profile <- detectHardware()
print(hw_profile)

# Configure parallel backend
backend <- configureParallelBackend(hw_profile)

# Use optimized thread count in analysis functions
fastq_dir <- "path/to/fastq/files"  # Replace with your actual path
fastq_files <- detectFastqFiles(fastq_dir, paired = TRUE)

# Run filtering with optimal thread count
filtered_files <- filterAndTrimReads(
  fastq_files, 
  output_dir = "filtered",
  truncLen = c(240, 200),
  maxEE = c(2, 2),
  multithread = hw_profile$optimal_threads  # Use detected optimal thread count
)

# Memory-aware processing for large datasets
error_memory <- estimateMemoryPerTask("error_learning", read_size = 1e9)
thread_count <- min(hw_profile$optimal_threads, 
                   floor(hw_profile$free_memory * 0.8 / error_memory))

# Run DADA2 pipeline with optimized parameters
results <- runDADA2Pipeline(
  filtered_files,
  min_length = 240,
  max_length = 260,
  multithread = thread_count  # Use memory-aware thread count
)
```

For detailed examples on hardware optimization, see the hardware vignette:

```r
vignette("hardware_optimization", package = "ampliSeqR")
```

The package leverages hardware acceleration through optimized parallel backends, DADA2's multithreading capabilities, and Rcpp integration for performance-critical operations.

## Advanced File Detection

ampliSeqR provides powerful file detection capabilities for various sequencing platforms and complex directory structures:

```r
library(ampliSeqR)
library(dplyr)  # For data manipulation

# Auto-detect platform with fuzzy matching (handles inconsistent naming)
fastq_files <- detectFastqFilesAdvanced(
  input_dir = "data/sequencing_run/",
  paired = TRUE,                     # Look for paired-end files
  recursive = TRUE,                  # Search in nested directories
  platform = "auto",                 # Auto-detect the sequencing platform
  fuzzy_threshold = 0.2,             # Allow fuzzy matching for filenames
  confidence_threshold = 0.7,        # Confidence threshold for matches
  extract_metadata = TRUE,           # Extract metadata from filenames
  verify_content = TRUE              # Verify actual FASTQ content
)

# Handle PacBio data with custom patterns
pacbio_files <- detectFastqFilesAdvanced(
  input_dir = "data/pacbio/",
  platform = "pacbio",
  confidence_threshold = 0.6  # Lower threshold for ambiguous naming
)

# Custom patterns for unusual naming conventions
custom_patterns <- list(
  forward_patterns = c("_forward_", "_f_", "_read1_"),
  reverse_patterns = c("_reverse_", "_r_", "_read2_"),
  sample_patterns = c(
    "(.+?)_forward_",
    "(.+?)_reverse_",
    "(.+?)_(f|r)_"
  )
)

# Use custom patterns
custom_files <- detectFastqFilesAdvanced(
  input_dir = "data/custom_naming/",
  paired = TRUE,
  known_patterns = custom_patterns
)

# Extract and use metadata from filenames
metadata_summary <- fastq_files %>%
  mutate(
    lane = sapply(metadata, function(x) ifelse(is.null(x$lane), NA, x$lane)),
    date = sapply(metadata, function(x) ifelse(is.null(x$date), NA, x$date)),
    flowcell = sapply(metadata, function(x) ifelse(is.null(x$flowcell), NA, x$flowcell))
  ) %>%
  select(sample_name, lane, date, flowcell, confidence)

# Filter to high-confidence matches only
high_confidence_pairs <- fastq_files %>% 
  filter(confidence > 0.9) %>%
  arrange(desc(confidence))

# Multi-platform study example
# Process each platform separately, then combine results
platforms <- c("illumina", "pacbio", "nanopore")
platform_results <- list()

for (platform_name in platforms) {
  # Detect files for this platform
  platform_files <- detectFastqFilesAdvanced(
    input_dir = paste0("data/", platform_name),
    platform = platform_name,
    recursive = TRUE
  )
  
  # Store results for each platform
  platform_results[[platform_name]] <- platform_files
  
  # Output some information about detected files
  cat(sprintf("Detected %d samples from %s platform\n", 
              nrow(platform_files), platform_name))
}

# Examine header formats from content verification
illumina_headers <- unique(platform_results[["illumina"]]$forward_header)
pacbio_headers <- unique(platform_results[["pacbio"]]$forward_header)
```

For more examples, see the included example script:

```r
file.edit(system.file("examples/advanced_file_detection.R", package = "ampliSeqR"))
```

### Integration with the Analysis Pipeline

The advanced file detection can be seamlessly integrated with the rest of the ampliSeqR pipeline:

```r
library(ampliSeqR)

# Use advanced file detection with a multi-batch study
fastq_files <- detectFastqFilesAdvanced(
  input_dir = "data/microbiome_study/",
  recursive = TRUE,  # Search in all subdirectories
  platform = "auto",  # Auto-detect the platform
  extract_metadata = TRUE
)

# Extract batch information from directory structure
fastq_files <- fastq_files %>%
  mutate(
    batch = basename(dirname(forward)),
    sequencing_date = sapply(metadata, function(x) ifelse(is.null(x$date), NA, x$date))
  )

# Process each batch separately with appropriate parameters
batches <- unique(fastq_files$batch)

results_list <- list()
for (batch_id in batches) {
  message("Processing batch: ", batch_id)
  
  # Filter to current batch
  batch_files <- fastq_files %>% filter(batch == batch_id)
  
  # Get optimal filtering parameters for this batch
  truncLen <- if(grepl("illumina", tolower(batch_id))) c(240, 200) else c(200, 150)
  
  # Filter and trim
  filtered_files <- filterAndTrimReads(
    batch_files,
    output_dir = file.path("filtered", batch_id),
    truncLen = truncLen,
    maxEE = c(2, 2)
  )
  
  # Run DADA2 pipeline
  batch_results <- runDADA2Pipeline(
    filtered_files,
    remove_chimeras = TRUE,
    multithread = TRUE
  )
  
  # Store results
  results_list[[batch_id]] <- batch_results
}

# Combine results from all batches for downstream analysis
combined_asv_table <- do.call(cbind, lapply(results_list, function(x) x$asv_table))
```

## Documentation

For detailed documentation on each function:

```r
?detectHardware           # Hardware detection and optimization
?configureParallelBackend # Parallel backend configuration
?estimateMemoryPerTask    # Memory requirement estimation
?detectFastqFiles         # Basic file detection and organization
?detectFastqFilesAdvanced # Advanced file detection with fuzzy matching
?filterAndTrimReads       # Read filtering and trimming
?runDADA2Pipeline         # DADA2 pipeline wrapper
?plotQualityProfiles      # Quality visualization
```

## License

This package is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use ampliSeqR in your research, please cite:

```
Developer, S. (2025). ampliSeqR: Accelerated 16S Amplicon Sequencing Analysis. R package version 0.1.0.
```
