# Advanced FASTQ File Detection Example
# This script demonstrates the use of the enhanced file detection capabilities
# in ampliSeqR for handling complex file naming conventions and directory structures.

library(ampliSeqR)

# Basic usage with auto-detected platform
# This will scan recursively and detect platform automatically
fastq_files <- detectFastqFilesAdvanced(
  input_dir = "data/",
  paired = TRUE,
  recursive = TRUE
)

# Print summary of detected files
print(fastq_files)

# Example with specific platform and lower confidence threshold
# Useful for ambiguous file naming or non-standard conventions
fastq_files_nanopore <- detectFastqFilesAdvanced(
  input_dir = "data/nanopore/",
  paired = TRUE,
  platform = "nanopore",
  confidence_threshold = 0.5,
  recursive = TRUE
)

# Example with custom patterns for unusual naming conventions
custom_patterns <- list(
  forward_patterns = c("_forward_", "_f_", "_read1_"),
  reverse_patterns = c("_reverse_", "_r_", "_read2_"),
  sample_patterns = c(
    "(.+?)_forward_",
    "(.+?)_reverse_",
    "(.+?)_(f|r)_"
  )
)

fastq_files_custom <- detectFastqFilesAdvanced(
  input_dir = "data/custom_naming/",
  paired = TRUE,
  known_patterns = custom_patterns,
  confidence_threshold = 0.6
)

# Example for multi-batch experiment with nested directories
# This will analyze the directory structure to identify batches
fastq_files_multi_batch <- detectFastqFilesAdvanced(
  input_dir = "data/experiment/",
  paired = TRUE,
  recursive = TRUE,
  extract_metadata = TRUE
)

# Example with content verification
# This will sample sequence headers to validate file content
fastq_files_verified <- detectFastqFilesAdvanced(
  input_dir = "data/",
  paired = TRUE,
  verify_content = TRUE, 
  verify_lines = 2000  # Check more lines for better validation
)

# Advanced filtering by confidence score
library(dplyr)

high_confidence_pairs <- fastq_files %>% 
  filter(confidence > 0.9) %>%
  arrange(desc(confidence))

# Extract metadata from detected files
metadata_summary <- fastq_files %>%
  mutate(
    lane = sapply(metadata, function(x) ifelse(is.null(x$lane), NA, x$lane)),
    date = sapply(metadata, function(x) ifelse(is.null(x$date), NA, x$date)),
    flowcell = sapply(metadata, function(x) ifelse(is.null(x$flowcell), NA, x$flowcell))
  ) %>%
  select(sample_name, lane, date, flowcell, confidence)

# Proceed with downstream analysis using the detected files
if (nrow(fastq_files) > 0) {
  # Example: Filter and trim the detected reads
  filtered_reads <- filterAndTrimReads(
    fastq_data = fastq_files,
    output_dir = "filtered/",
    truncLen = c(240, 200),
    maxEE = c(2, 2)
  )
}