






plot_flowset_histograms <- function(
    fs,
    gate = NULL,
    channel = "FSC-A",
    transformation = c("none", "asinh"),
    cofactor = 150,
    bins = 200,
    xlab = NULL,
    ylab = "Density"
) {
  
  transformation <- match.arg(transformation)
  
  # Validate flowSet
  if (!inherits(fs, "flowSet")) {
    stop("'fs' must be a flowSet object.", call. = FALSE)
  }

  # Check that the requested channel is present in all samples
  has_channel <- fsApply(fs, function(ff) {sum(channel == colnames(ff))})

  if (!all(has_channel == 1)) {
    missing_samples <- sampleNames(fs)[has_channel != 1]
    
    stop(
      sprintf(
        "Channel '%s' is missing from: %s",
        channel,
        paste(missing_samples, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Validate transformation-specific arguments
  if (transformation == "asinh" &&
      (!is.numeric(cofactor) || length(cofactor) != 1 ||
       !is.finite(cofactor) || cofactor <= 0)) {
    stop(
      "'cofactor' must be a single positive finite number when ",
      "transformation = 'asinh'.",
      call. = FALSE
    )
  }

  # Validate number of bins
  if (!is.numeric(bins) || length(bins) != 1 ||
      !is.finite(bins) || bins <= 0) {
    stop("'bins' must be a positive number.", call. = FALSE)
  }

  # Apply gate
  if (!is.null(gate)) {
    fs <- flowCore::Subset(
      fs,
      gate
    )
  }

  # Default x-axis label
  if (is.null(xlab)) {
    xlab <- channel
  }

  # Generate one plot per sample
  plots <- lapply(seq_along(fs), function(i) {
    
    ff <- fs[[i]]
    signal <- flowCore::exprs(ff)[, channel]

    if (transformation == "asinh") {
      signal <- asinh(signal / cofactor)
    }

    df <- data.frame(signal = signal)

    ggplot2::ggplot( df, aes(x = signal)) +
       ggplot2::geom_histogram(bins = bins, aes(y = after_stat(density))) +
       ggplot2::labs(title = sampleNames(fs)[i], x = xlab,y = ylab) +
      cowplot::theme_cowplot()

  })

  names(plots) <- sampleNames(fs)

  plots
}

