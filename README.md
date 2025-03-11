# ampliSeqR

An R package for efficient 16S amplicon sequencing analysis with hardware acceleration.

## Overview

ampliSeqR provides tools for analyzing 16S rRNA gene amplicon sequencing data, with a focus on performance and hardware acceleration. This package wraps and extends the functionality of DADA2 to simplify and optimize the analysis workflow.

## Features

- Automatic detection and organization of FASTQ files
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
fastq_files <- detectFastqFiles(fastq_dir, paired = TRUE)
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

The package leverages hardware acceleration through DADA2's multithreading capabilities and Rcpp integration for performance-critical operations. Use the `multithread` parameter when available to take advantage of multiple CPU cores.

## Documentation

For detailed documentation on each function:

```r
?detectFastqFiles
?filterAndTrimReads
?runDADA2Pipeline
?plotQualityProfiles
```

## License

This package is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use ampliSeqR in your research, please cite:

```
Developer, S. (2025). ampliSeqR: Accelerated 16S Amplicon Sequencing Analysis. R package version 0.1.0.
```