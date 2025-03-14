% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plotQualityProfiles}
\alias{plotQualityProfiles}
\title{Plot quality profiles}
\usage{
plotQualityProfiles(
  quality_data,
  aggregate = TRUE,
  n_samples = 4,
  read_type = "both"
)
}
\arguments{
\item{quality_data}{Quality data as produced by calculateQualityStats()}

\item{aggregate}{Logical indicating whether to aggregate samples (default: TRUE)}

\item{n_samples}{Integer indicating number of samples to plot if not aggregating (default: 4)}

\item{read_type}{Character indicating which read type to plot ("forward", "reverse", or "both") (default: "both")}
}
\value{
A ggplot object with quality profile(s)
}
\description{
This function plots quality profiles for FastQ files using ggplot2.
}
\examples{
\dontrun{
fastq_files <- detectFastqFiles("path/to/fastq/dir")
quality_data <- calculateQualityStats(fastq_files)
plotQualityProfiles(quality_data)
}
}
