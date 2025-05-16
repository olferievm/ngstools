#'Name: plot_2D_PCA
#'
#'Purpose:
#  Generates a 2D PCA (or MDS) scatter plot using the top principal components or dimensions, color-coded by group, with flexible options for aesthetics, labeling, axis reversal, and subsetting.

#Parameters:
#  Argument	Description
#' @param x - PCA or MDS object (prcomp or custom MDS with $eigen.vectors, $eigen.values)
#' @param group - Grouping variable (factor or numeric)
#' @param ind - Labels for individual points
#' @param pc - Integer vector of length 2 specifying which components to plot
#' @param col - Custom colors for each group (optional)
#' @param alpha - Transparency of points (default 1)
#' @param xlim 	 Optional X axis limits; computed automatically if NULL
#' @param ylim 	 Optional Y axis limits; computed automatically if NULL
#' @param subgroup 	 Logical or integer index to subset the samples to plot
#' @param main 	 Plot title (currently unused)
#' @param col.text - Color of text labels
#' @param text.size - Font size for point labels
#' @param point.size - Size of points on the plot
#' @param segment.color - Color of text label leader lines
#' @param revX - Logical to reverse y axes
#' @param revY - Logical to reverse y axes
#' @param show.labels - Whether to display point labels
#' @param legend.title - Title for the legend



plot_2D_PCA <- function(x,
                        group,
                        ind,
                        pc = c(1, 2),
                        col = NULL,
                        alpha = 1,
                        xlim = NULL,
                        ylim = NULL,
                        subgroup = NULL,
                        main = "Principal components",
                        col.text = "black",
                        text.size = 2,
                        point.size = 5,
                        segment.color = "#555555",
                        revX = FALSE,
                        revY = FALSE,
                        show.labels = TRUE,
                        legend.title = "Group") {
  # Check input object
  if (!class(x)[1] %in% c("MDS", "prcomp")) {
    stop("Use MDS or prcomp PCA object")
  }
  
  # Load required packages
  require(ggplot2)
  require(RColorBrewer)
  require(ggrepel)
  require(cowplot)
  
  # Check valid PC components
  if (!is.numeric(pc) || length(pc) != 2) {
    stop("Specify the selected components with an integer vector of length 2")
  }
  
  if (length(group) != length(ind)) {
    stop("Length of group and ind differs")
  }
  
  # Prepare color palette based on factor levels
  if (is.character(group)) group <- factor(group)
  
  if (is.factor(group)) {
    Ngroup <- length(levels(group))
    
    if (!is.null(col)) {
      if (length(col) > Ngroup) {
        col <- col[1:Ngroup]
      }
    }
    
    if (is.null(col) || length(col) != Ngroup) {
      col <- if (Ngroup <= 8) {
        brewer.pal(n = Ngroup, name = "Dark2")
      } else {
        sample(colorRampPalette(brewer.pal(n = 12, name = "Paired"))(Ngroup))
      }
    }
  }
  
  # Internal helper function: scale limits to nearest quarter (0.25) range
  quaterScale <- function(x) {
    logll <- log10(abs(x))
    dr <- floor(logll)
    x <- x / (10^dr)
    ds <- floor(x)
    dp <- x - trunc(x)
    d <- ifelse(dp < 0.25, ds + 0.25,
                ifelse(dp < 0.5, ds + 0.5,
                       ifelse(dp < 0.75, ds + 0.75,
                              ds + 1)))
    d * 10^dr
  }
  
  # Calculate optimal axis limits with padding
  opt_limits <- function(x, f = 0.05) {
    ll <- min(x)
    ll <- ifelse(ll >= 0, ll - (ll * f), ll + (ll * f))
    if (ll != 0) {
      sign_ll <- ifelse(ll < 0, -1, 1)
      ll <- quaterScale(abs(ll)) * sign_ll
    }
    
    ul <- max(x)
    ul <- ifelse(ul >= 0, ul + (ul * f), ul - (ul * f))
    if (ul != 0) {
      sign_ul <- ifelse(ul < 0, -1, 1)
      ul <- quaterScale(abs(ul)) * sign_ul
    }
    
    return(c(Low.Limit = ll, Up.Limit = ul))
  }
  
  # Extract PCA values (assumes prcomp or MDS object structure)
  pc1 <- x$eigen.vectors[, pc[1]]
  pc2 <- x$eigen.vectors[, pc[2]]
  
  v1 <- round(x$eigen.values[pc[1]] / sum(x$eigen.values) * 100, 0)
  v2 <- round(x$eigen.values[pc[2]] / sum(x$eigen.values) * 100, 0)
  
  xLab <- sprintf("Dim%02d %s%%", pc[1], v1)
  yLab <- sprintf("Dim%02d %s%%", pc[2], v2)
  
  # Build data frame for plotting
  plot_df <- data.frame(
    PC1 = pc1,
    PC2 = pc2,
    Group = group,
    ind = ind,
    stringsAsFactors = FALSE
  )
  
  # Apply subgrouping if requested
  if (!is.null(subgroup)) {
    plot_df <- droplevels(plot_df[subgroup, ])
  }
  
  # Compute axis limits if not provided
  if (is.null(xlim)) xlim <- opt_limits(plot_df$PC1)
  if (is.null(ylim)) ylim <- opt_limits(plot_df$PC2)
  
  # Start building ggplot
  g <- ggplot(plot_df, aes(x = PC1, y = PC2, group = Group, color = Group)) +
    geom_point(size = point.size, alpha = alpha)
  
  # Apply custom or default color scale
  if (is.factor(group)) {
    g <- g + scale_color_manual(labels = levels(plot_df$Group), values = col)
  } else if (is.numeric(group)) {
    g <- g + scale_color_viridis_c()
  }
  
  # Add labels, limits, and legend
  g <- g + labs(x = xLab, y = yLab, color = legend.title)
  
  if (revX) {
    g <- g + xlim(rev(xlim))
  } else {
    g <- g + xlim(xlim)
  }
  
  if (revY) {
    g <- g + ylim(rev(ylim))
  } else {
    g <- g + ylim(ylim)
  }
  
  # Add text labels if requested
  if (show.labels) {
    g <- g + geom_text_repel(
      data = plot_df,
      aes(label = ind),
      color = col.text,
      size = text.size,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"),
      segment.color = segment.color,
      show.legend = NA
    )
  }
  
  return(g)
}


