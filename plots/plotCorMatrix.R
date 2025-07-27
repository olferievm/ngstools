#' Plot Annotated Correlation Matrix as Heatmap
#'
#' @description
#' This function plots a correlation matrix as a heatmap using ggplot2,
#' with optional overlaid labels for p-values (or other associated statistics).
#'
#' @param x A numeric matrix of correlation values (e.g., Pearson/Spearman r).
#' @param p A numeric matrix of corresponding p-values (same dimensions as x).
#' @param col_palette Vector of colors for gradient fill. Defaults to a Blue-White-Red palette.
#' @param text.size Size of text labels for p-values.
#' @param axis.text.size Size of axis text (row/column names).
#' @param text.color Color of overlaid text labels.
#' @param label.precision Decimal places to round p-values for display.
#' @param label.threshold Threshold for displaying p-values. Values above are omitted.
#' @param reorder.rows Logical. Cluster and reorder rows using hierarchical clustering.
#' @param reorder.cols Logical. Cluster and reorder columns using hierarchical clustering.
#' @param xmin Minimum value for fill scale (default = -1).
#' @param xmid Midpoint for fill scale (default = 0).
#' @param xmax Maximum value for fill scale (default = 1).
#' @param main Main title for the plot.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#'
#' @return A ggplot2 heatmap object.
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- matrix(runif(9, -1, 1), 3, dimnames = list(c("A", "B", "C"), c("X", "Y", "Z")))
#' p <- matrix(runif(9, 0, 0.1), 3, dimnames = list(c("A", "B", "C"), c("X", "Y", "Z")))
#' plotCorMatrix(x, p)
plotCorMatrix <- function(x, p, 
                          col_palette = NULL,
                          text.size = 3, 
                          axis.text.size = 10, 
                          text.color = 'green', 
                          label.precision = 3, 
                          label.threshold = 0.05, 
                          reorder.rows = TRUE,
                          reorder.cols = TRUE,
                          xmin = -1, xmid = 0, xmax = 1,
                          main = "", xlab = "", ylab = "") {
  
  # --- Input checks ---
  if (!is.matrix(x) || !is.numeric(x)) stop("x must be a numeric matrix.")
  if (!is.matrix(p) || !is.numeric(p)) stop("p must be a numeric matrix.")
  if (!all(dim(x) == dim(p))) stop("x and p must have the same dimensions.")
  
  # --- Default palette if none provided ---
  if (is.null(col_palette)) {
    col_palette <- c("#0803A1", "#1514A6", "#225EA8", "#74A9CF", "#BDC9E1", "#F1EEF6",
                     "#FDCC8A", "#FC8D59", "#D7301F", "#C11C38", "#b01919")
  }
  
  pcut <- 10^(-label.precision)
  # --- Prepare data ---
  df_x <- reshape2::melt(x, varnames = c("cluster", "group"), value.name = "value")
  df_p <- reshape2::melt(p, varnames = c("cluster", "group"), value.name = "label")
  df <- dplyr::left_join(df_x, df_p, by = c("cluster", "group")) %>%
    dplyr::mutate(
      labelm = case_when(
        is.na(label) ~ "",
        label < label.threshold & label >= pcut ~ as.character(round(label, label.precision)),
        label < label.threshold & label < pcut ~ sprintf("< %s", pcut),
        label > label.threshold ~ "",
        TRUE ~ ""
    ))
  
  # --- Reorder rows/columns if needed ---
  if (reorder.rows && nrow(x) > 1) {
    h <- hclust(dist(x), method = "ward.D2")
    df$cluster <- factor(df$cluster, levels = rownames(x)[h$order])
  }
  if (reorder.cols && ncol(x) > 1) {
    h <- hclust(dist(t(x)), method = "ward.D2")
    df$group <- factor(df$group, levels = colnames(x)[h$order])
  }
  
  # --- Heatmap plot ---
  g <- ggplot(df, aes(x = group, y = cluster, fill = value, label = labelm)) +
    geom_tile() +
    geom_text(data = dplyr::filter(df, !is.na(labelm), labelm != ""), 
              size = text.size, color = text.color) +
    scale_fill_gradientn(
      colours = col_palette,
      limits = c(xmin, xmax),
      rescaler = ~ scales::rescale_mid(.x, mid = xmid)
    ) +
    labs(title = main, x = xlab, y = ylab) +
    theme_minimal() +
    theme(axis.text = element_text(size = axis.text.size))
  
  return(g)
}
