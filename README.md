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

# Install from GitHub
devtools::install_github("yourusername/ampliSeqR")
```

## Quick Start

```r
library(ampliSeqR)

# Detect FASTQ files
fastq_files <- detectFastqFiles("path/to/fastq/dir", paired = TRUE)

# Analyze quality
quality_data <- calculateQualityStats(fastq_files)
plotQualityProfiles(quality_data)

# Filter and trim reads
filtered_files <- filterAndTrimReads(
  fastq_files, 
  "path/to/filtered/dir",
  truncLen = c(240, 200),
  maxEE = c(2, 2)
)

# Run DADA2 pipeline
results <- runDADA2Pipeline(filtered_files)

# Visualize results
plotASVAbundances(results)
```

## Documentation

See the package vignettes for more detailed examples and workflows:

```r
browseVignettes("ampliSeqR")
```

## License

This package is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use ampliSeqR in your research, please cite:

```
Citation information will be available here.
```