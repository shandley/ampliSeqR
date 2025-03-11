#' Plot quality profiles
#'
#' This function plots quality profiles for FastQ files using ggplot2.
#'
#' @param quality_data Quality data as produced by calculateQualityStats()
#' @param aggregate Logical indicating whether to aggregate samples (default: TRUE)
#' @param n_samples Integer indicating number of samples to plot if not aggregating (default: 4)
#' @param read_type Character indicating which read type to plot ("forward", "reverse", or "both") (default: "both")
#'
#' @return A ggplot object with quality profile(s)
#'
#' @export
#' @importFrom dplyr %>% filter group_by summarize
#' @importFrom tidyr unnest
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon theme_bw facet_wrap labs scale_color_brewer
#' @importFrom scales percent
#'
#' @examples
#' \dontrun{
#' fastq_files <- detectFastqFiles("path/to/fastq/dir")
#' quality_data <- calculateQualityStats(fastq_files)
#' plotQualityProfiles(quality_data)
#' }
plotQualityProfiles <- function(quality_data,
                               aggregate = TRUE,
                               n_samples = 4,
                               read_type = "both") {
  
  # Input validation
  if (is.null(quality_data$forward)) {
    stop("Forward quality data not found")
  }
  
  if (read_type == "both" && is.null(quality_data$reverse)) {
    warning("Reverse quality data not found, showing only forward")
    read_type <- "forward"
  }
  
  # Process forward reads
  forward_data <- quality_data$forward %>%
    mutate(read_type = "Forward")
  
  # Process reverse reads if available and requested
  if (read_type %in% c("reverse", "both") && !is.null(quality_data$reverse)) {
    reverse_data <- quality_data$reverse %>%
      mutate(read_type = "Reverse")
    
    if (read_type == "both") {
      plot_data <- bind_rows(forward_data, reverse_data)
    } else {
      plot_data <- reverse_data
    }
  } else {
    plot_data <- forward_data
  }
  
  # Extract position-specific quality data
  plot_data <- plot_data %>%
    mutate(
      quality_data = map(quality, function(x) {
        # Check if quality is null (failed quality check)
        if (is.null(x)) return(NULL)
        
        # Extract per-cycle quality
        per_cycle <- x[["perCycle"]]$quality
        if (is.null(per_cycle)) return(NULL)
        
        # Convert to tibble for easier manipulation
        tibble(
          Cycle = per_cycle$Cycle,
          Score = per_cycle$Score
        )
      })
    ) %>%
    filter(!map_lgl(quality_data, is.null))
  
  # Sample n_samples if not aggregating
  if (!aggregate && nrow(plot_data) > n_samples) {
    set.seed(42)  # For reproducibility
    sample_idx <- sample(1:nrow(plot_data), n_samples)
    plot_data <- plot_data[sample_idx, ]
  }
  
  # Prepare data for plotting
  if (aggregate) {
    # Aggregate across samples
    plot_quality <- plot_data %>%
      select(sample_name, read_type, quality_data) %>%
      unnest(quality_data) %>%
      group_by(read_type, Cycle) %>%
      summarize(
        mean_quality = mean(Score),
        min_quality = min(Score),
        q25 = quantile(Score, 0.25),
        median_quality = median(Score),
        q75 = quantile(Score, 0.75),
        max_quality = max(Score),
        .groups = "drop"
      )
    
    # Create plot
    p <- ggplot(plot_quality, aes(x = Cycle, y = mean_quality, color = read_type)) +
      geom_ribbon(aes(ymin = q25, ymax = q75, fill = read_type), alpha = 0.2, color = NA) +
      geom_line(size = 1) +
      labs(
        x = "Cycle (bp position)",
        y = "Quality Score",
        title = "Quality Profile",
        subtitle = paste("Aggregated across", nrow(plot_data), "samples")
      ) +
      scale_color_brewer(palette = "Set1") +
      scale_fill_brewer(palette = "Set1") +
      theme_bw() +
      theme(legend.title = element_blank())
    
  } else {
    # Plot individual samples
    plot_quality <- plot_data %>%
      select(sample_name, read_type, quality_data) %>%
      unnest(quality_data)
    
    # Create plot
    p <- ggplot(plot_quality, aes(x = Cycle, y = Score, color = read_type)) +
      geom_line(alpha = 0.8) +
      facet_wrap(~sample_name) +
      labs(
        x = "Cycle (bp position)",
        y = "Quality Score",
        title = "Quality Profiles by Sample"
      ) +
      scale_color_brewer(palette = "Set1") +
      theme_bw() +
      theme(legend.title = element_blank())
  }
  
  return(p)
}

#' Plot read length distribution
#'
#' This function plots the distribution of read lengths.
#'
#' @param length_data Length distribution data as produced by analyzeReadLengths()
#' @param aggregate Logical indicating whether to aggregate samples (default: TRUE)
#' @param n_samples Integer indicating number of samples to plot if not aggregating (default: 4)
#'
#' @return A ggplot object with length distribution(s)
#'
#' @export
#' @importFrom dplyr %>% filter group_by summarize
#' @importFrom tidyr unnest
#' @importFrom ggplot2 ggplot aes geom_histogram theme_bw facet_wrap labs
#'
#' @examples
#' \dontrun{
#' fastq_files <- detectFastqFiles("path/to/fastq/dir")
#' length_data <- analyzeReadLengths(fastq_files)
#' plotReadLengths(length_data)
#' }
plotReadLengths <- function(length_data,
                           aggregate = TRUE,
                           n_samples = 4) {
  
  # Input validation
  if (nrow(length_data) == 0) {
    stop("No length data provided")
  }
  
  # Sample n_samples if not aggregating
  if (!aggregate && nrow(length_data) > n_samples) {
    set.seed(42)  # For reproducibility
    sample_idx <- sample(1:nrow(length_data), n_samples)
    plot_data <- length_data[sample_idx, ]
  } else {
    plot_data <- length_data
  }
  
  # Unnest the length distributions
  plot_lengths <- plot_data %>%
    select(sample_name, read_type, length_dist) %>%
    unnest(length_dist)
  
  # Create plot
  if (aggregate) {
    # Aggregate across samples
    agg_lengths <- plot_lengths %>%
      group_by(read_type, length) %>%
      summarize(
        total_count = sum(count),
        .groups = "drop"
      ) %>%
      group_by(read_type) %>%
      mutate(proportion = total_count / sum(total_count))
    
    p <- ggplot(agg_lengths, aes(x = length, y = proportion, fill = read_type)) +
      geom_col(position = "dodge") +
      labs(
        x = "Read Length (bp)",
        y = "Proportion of Reads",
        title = "Read Length Distribution",
        subtitle = paste("Aggregated across", nrow(plot_data), "samples")
      ) +
      scale_fill_brewer(palette = "Set1") +
      theme_bw() +
      theme(legend.title = element_blank())
    
  } else {
    # Plot individual samples
    p <- ggplot(plot_lengths, aes(x = length, y = proportion, fill = read_type)) +
      geom_col(position = "dodge") +
      facet_wrap(~sample_name) +
      labs(
        x = "Read Length (bp)",
        y = "Proportion of Reads",
        title = "Read Length Distribution by Sample"
      ) +
      scale_fill_brewer(palette = "Set1") +
      theme_bw() +
      theme(legend.title = element_blank())
  }
  
  return(p)
}

#' Plot ASV abundance distribution
#'
#' This function plots the distribution of ASV abundances.
#'
#' @param asv_results ASV results as produced by generateASVTable() or runDADA2Pipeline()
#' @param top_n Integer indicating number of top ASVs to highlight (default: 10)
#' @param log_scale Logical indicating whether to use log scale for abundance (default: TRUE)
#'
#' @return A ggplot object with ASV abundance distribution
#'
#' @export
#' @importFrom dplyr %>% arrange desc slice mutate
#' @importFrom ggplot2 ggplot aes geom_point theme_bw scale_y_log10 labs scale_color_manual scale_size_manual element_text
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' \dontrun{
#' results <- runDADA2Pipeline(filtered_files)
#' plotASVAbundances(results)
#' }
plotASVAbundances <- function(asv_results, top_n = 10, log_scale = TRUE) {
  
  # Input validation
  if (is.null(asv_results$asv_stats)) {
    stop("ASV statistics not found in results")
  }
  
  # Prepare data for plotting
  asv_stats <- asv_results$asv_stats
  
  # Identify top ASVs by abundance
  top_asvs <- asv_stats %>%
    arrange(desc(abundance)) %>%
    slice(1:min(top_n, nrow(asv_stats))) %>%
    pull(asv_id)
  
  # Create plotting data with highlighted top ASVs
  plot_data <- asv_stats %>%
    mutate(
      highlight = asv_id %in% top_asvs,
      label = if_else(highlight, asv_id, NA_character_)
    )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = prevalence, y = abundance, color = highlight, label = label)) +
    geom_point(aes(size = highlight), alpha = 0.7) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50")) +
    scale_size_manual(values = c("TRUE" = 3, "FALSE" = 1)) +
    ggrepel::geom_text_repel(
      data = filter(plot_data, highlight),
      box.padding = 0.5,
      max.overlaps = Inf,
      size = 3
    ) +
    labs(
      x = "Prevalence (number of samples)",
      y = "Abundance (number of reads)",
      title = "ASV Abundance Distribution",
      subtitle = paste("Top", top_n, "ASVs highlighted")
    ) +
    theme_bw() +
    theme(legend.position = "none")
  
  # Apply log scale if requested
  if (log_scale) {
    p <- p + scale_y_log10()
  }
  
  return(p)
}

