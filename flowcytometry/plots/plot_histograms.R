#' Plot flow cytometry channel histograms
#'
#' Generates density histograms for selected channels in a \code{flowFrame}.
#' An optional gate can be applied before plotting, and fluorescence values
#' can optionally be transformed using an inverse hyperbolic sine
#' transformation.
#'
#' If \code{channels = NULL}, all channels in the supplied \code{flowFrame}
#' are used except those matching \code{exclude_channels}.
#'
#' @param ff A \code{flowFrame} containing flow cytometry data.
#'
#' @param gate Optional \code{filter} object, such as a \code{polygonGate},
#'   used to subset events before plotting. If \code{NULL}, all events are
#'   retained.
#'
#' @param channels Character vector specifying channels to plot. If
#'   \code{NULL}, channels are selected automatically from \code{ff} after
#'   removing channels matching \code{exclude_channels}.
#'
#' @param exclude_channels Character vector of channel names or patterns to
#'   exclude when \code{channels = NULL}. Defaults to
#'   \code{c("Autofluo", "Time", "FSC-A", "SSC-A")}. Set to \code{NULL}
#'   to disable exclusion.
#'
#' @param transform Character string specifying the transformation applied
#'   before plotting. One of \code{"asinh"} or \code{"none"}.
#'
#' @param cofactor Numeric cofactor used for the \code{"asinh"}
#'   transformation. Values are transformed as
#'   \code{asinh(x / cofactor)}. Defaults to 150.
#'
#' @param bins Integer specifying the number of histogram bins. Defaults
#'   to 200.
#'
#' @param xlab Character string specifying the x-axis label.
#'
#' @param ylab Character string specifying the y-axis label.
#'
#' @return A named list of \code{ggplot} objects, with one histogram for
#'   each selected channel.
#'
#' @importFrom flowCore Subset exprs
#' @importFrom ggplot2 ggplot aes geom_histogram labs after_stat
#' @importFrom cowplot theme_cowplot
#' @importFrom stringr str_subset
#'
#' @examples
#' \dontrun{
#' p <- plot_histograms(
#'     fs_cells_unmix[[1]],
#'     gate = cell_gate,
#'     transform = "asinh",
#'     cofactor = 150
#' )
#'
#' p[["CD3_BUV805"]]
#' }
#'
#' @export


plot_histograms <- function(
    ff,
    gate = NULL,
    channels = NULL,
    exclude_channels = c("Autofluo", "Time", "FSC-A", "SSC-A"),
    transform = c("asinh", "none"),
    cofactor = 150,
    bins = 200,
    xlab = "Unmixed signal",
    ylab = "Density"
) {
  
  transform <- match.arg(transform)
  
  # Apply gate
  if (!is.null(gate)) {
    ff <- flowCore::Subset(ff, gate)
  }
  
  # Determine available channels
  available_channels <- colnames(ff)
  
  # Select channels automatically
  if (is.null(channels)) {
    
    channels <- available_channels
    
    if (!is.null(exclude_channels)) {
      pattern <- paste(exclude_channels, collapse = "|")
      
      channels <- stringr::str_subset(
        channels,
        pattern = pattern,
        negate = TRUE
      )
    }
    
  } else {
    
    # Check explicitly requested channels
    missing_channels <- setdiff(channels, available_channels)
    
    if (length(missing_channels) > 0) {
      stop(
        "Channels not found in flowFrame: ",
        paste(missing_channels, collapse = ", ")
      )
    }
  }
  
  # Extract expression matrix once
  x <- flowCore::exprs(ff)
  
  plots <- lapply(channels, function(ch) {
    
    signal <- x[, ch]
    
    if (transform == "asinh") {
      signal <- asinh(signal / cofactor)
    }
    
    df <- data.frame(signal = signal)
    
    ggplot2::ggplot(df, ggplot2::aes(x = signal)) +
      ggplot2::geom_histogram(
        bins = bins,
        ggplot2::aes(y = ggplot2::after_stat(density))
      ) +
      ggplot2::labs(
        title = ch,
        x = xlab,
        y = ylab
      ) +
      cowplot::theme_cowplot()
  })
  
  names(plots) <- channels
  
  plots
}