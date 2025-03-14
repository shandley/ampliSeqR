% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dada2_wrapper.R
\name{generateASVTable}
\alias{generateASVTable}
\title{Generate ASV table from DADA2 results}
\usage{
generateASVTable(
  dada_results,
  filtered_files,
  min_length = 0,
  max_length = Inf,
  remove_chimeras = TRUE,
  method = "consensus",
  verbose = TRUE,
  paired = TRUE
)
}
\arguments{
\item{dada_results}{List with dada objects as produced by inferSampleComposition()}

\item{filtered_files}{Tibble with filtered file information}

\item{min_length}{Minimum acceptable sequence length (default: 0)}

\item{max_length}{Maximum acceptable sequence length (default: Inf)}

\item{remove_chimeras}{Logical indicating whether to remove chimeras (default: TRUE)}

\item{method}{Character string specifying chimera removal method (default: "consensus")}

\item{verbose}{Logical indicating whether to print verbose output (default: TRUE)}

\item{paired}{Logical indicating whether data is paired-end (default: TRUE)}
}
\value{
A list containing the ASV table, ASV sequences, and summary statistics
}
\description{
This function generates an ASV (Amplicon Sequence Variant) table from the
results of DADA2 processing.
}
\examples{
\dontrun{
fastq_files <- detectFastqFiles("path/to/fastq/dir")
filtered_files <- filterAndTrimReads(fastq_files, "path/to/filtered/dir")
error_models <- learnErrorRates(filtered_files)
dada_objects <- inferSampleComposition(filtered_files, error_models)
asv_results <- generateASVTable(dada_objects, filtered_files)
}
}
