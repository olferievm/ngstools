#' Plot Biological Coefficient of Variation Using Plotly
#'
#' This function creates an interactive BCV (Biological Coefficient of Variation) plot
#' from a DGEList object using Plotly, showing tagwise, trended, and common dispersion estimates.
#'
#' @param y A `DGEList` object from the `edgeR` package containing dispersion estimates and gene annotation.
#' @param xlab A character string for the x-axis label. Default is "Average log CPM".
#' @param ylab A character string for the y-axis label. Default is "Biological coefficient of variation".
#' 
#' @return A plotly object displaying the BCV plot.
#' @import plotly
#' @export
#'
#' @examples
#' # Assuming 'dge' is a properly processed DGEList with dispersion estimates:
#' plotlyBCV(dge)

plotlyBCV <- function(y,
                      xlab = "Average log CPM",
                      ylab = "Biological coefficient of variation") {
  
  require(plotly, quietly = TRUE, warn.conflicts = FALSE)
  
  # Validate input
  if (!is(y, "DGEList")) 
    stop("y must be a DGEList.")
  
  # Use precomputed average log CPM or compute it from counts
  A <- y$AveLogCPM
  if (is.null(A)) 
    A <- aveLogCPM(y$counts, offset = getOffset(y))
  
  # Get dispersion values
  disp <- getDispersion(y)
  if (is.null(disp)) 
    stop("No dispersions to plot.")
  
  # If common dispersion, expand to vector
  if (attr(disp, "type") == "common") {
    disp <- rep(disp, length = length(A))
  }
  
  # Build data frame for plotting
  df <- data.frame(
    gene = y$genes$gene_name,
    chr = y$genes$chr,
    type = y$genes$gene_type,
    AverageCPM = A,
    Dispersion = sqrt(disp)  # BCV is the square root of dispersion
  )
  
  # Add optional dispersion estimates if available
  if (!is.null(y$tagwise.dispersion)) 
    df$Tagwise <- sqrt(y$tagwise.dispersion)
  if (!is.null(y$trended.dispersion)) 
    df$Trended <- sqrt(y$trended.dispersion)
  if (!is.null(y$common.dispersion)) 
    df$Common <- sqrt(y$common.dispersion)
  
  # Create interactive plot
  p <- plot_ly(df, x = ~AverageCPM)
  
  # Add tagwise (default) dispersion markers
  p <- p %>% add_markers(
    y = ~Dispersion,
    name = "Dispersion",
    type = 'scatter',
    mode = "markers",
    marker = list(color = "black"),
    text = ~paste(gene, '<br>', chr, '<br>', type)
  )
  
  # Add optional common dispersion line
  if (!is.null(y$common.dispersion)) {
    p <- p %>% add_lines(
      y = ~Common,
      name = 'Common',
      mode = 'lines',
      line = list(color = 'blue')
    )
  }
  
  # Add optional trended dispersion line
  if (!is.null(y$trended.dispersion)) {
    p <- p %>% add_lines(
      y = ~Trended,
      name = 'Trended',
      mode = 'lines',
      line = list(color = 'red')
    )
  }
  
  # Final layout adjustments
  p <- p %>% layout(
    title = "Dispersion",
    xaxis = list(title = xlab),
    yaxis = list(title = ylab, showgrid = FALSE)
  )
  
  return(p)
}
