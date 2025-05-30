#' Volcano Plot for Gene Expression Results
#'
#' This function generates a volcano plot from a `TopTags` object (from edgeR)
#' or a data frame containing gene expression results with `logFC`, `FDR`, and
#' a gene label column.
#'
#' @param top A `TopTags` object or a data frame with columns `logFC`, `FDR`, and `label_column`.
#' @param logFC Numeric. Threshold for absolute log2 fold change. Default: 1.
#' @param pval Numeric. FDR cutoff for significance (between 0 and 1). Default: 0.05.
#' @param label Logical. Whether to label significant genes. Default: TRUE.
#' @param main Character. Plot title. Default: "Volcano plot".
#' @param xlim "auto" or numeric vector of length 2 for x-axis limits. Default: "auto".
#' @param ylim "auto" or numeric vector of length 2 for y-axis limits. Default: "auto".
#' @param max.overlaps Integer. Max overlapping labels (ggrepel). Default: 10.
#' @param label_column Character. Column name to use for text labels. Default: "gene_name".
#' @param color_scale Character vector of length 2. Colors for c("Not Signif.", significant_label). Default: c("black", "orange").
#' @param label_size Numeric. Text size for gene labels. Default: 3.
#'
#' @return A `ggplot2` object.
#'
#' @import ggplot2
#' @import ggrepel
#' @import dplyr
#' @export
volcano_plot <- function(top,
                         logFC = 1,
                         pval = 0.05,
                         label = TRUE,
                         main = "Volcano plot",
                         xlim = "auto",
                         ylim = "auto",
                         max.overlaps = 10,
                         label_column = "gene_name",
                         color_scale = c("black", "orange"),
                         label_size = 3) {
  
  # Load required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Package 'ggrepel' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  # Convert from TopTags if needed
  if (inherits(top, "TopTags")) top <- top$table
  if (!is.data.frame(top)) stop("'top' must be a data.frame or a TopTags object.")
  
  # Check required columns
  required_cols <- c("logFC", "FDR", label_column)
  missing_cols <- setdiff(required_cols, colnames(top))
  if (length(missing_cols) > 0)
    stop("Missing column(s): ", paste(missing_cols, collapse = ", "))
  
  # Check input validity
  if (!(logFC > 0)) stop("logFC must be > 0.")
  if (pval <= 0 || pval >= 1) stop("pval must be between 0 and 1.")
  if (!is.character(color_scale) || length(color_scale) != 2)
    stop("color_scale must be a character vector of length 2.")
  
  # Axis limits
  if (is.character(xlim) && xlim == "auto") {
    xmax <- ceiling(max(abs(top$logFC), na.rm = TRUE)) + 1
    xlim <- c(-xmax, xmax)
  }
  if (is.character(ylim) && ylim == "auto") {
    ymax <- ceiling(max(-log10(top$FDR), na.rm = TRUE)) + 1
    ylim <- c(0, ymax)
  }
  
  # Define dynamic significance label
  signif_label <- sprintf("FDR < %.2g", pval)
  
  # Annotate data with significance
  top$Threshold <- ifelse(abs(top$logFC) > logFC & top$FDR < pval,
                          signif_label, "Not Signif.")
  
  # Subset significant points for labeling
  topnames <- dplyr::filter(top, Threshold == signif_label)
  
  # Create plot
  g <- ggplot(top, aes(x = logFC, y = -log10(FDR), color = Threshold)) +
    geom_point(alpha = 0.5, size = 1.75) +
    xlab("log2 fold change") +
    ylab("-log10 FDR") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom") +
    xlim(xlim) +
    ylim(ylim) +
    ggtitle(main) +
    scale_color_manual(values = setNames(color_scale, c("Not Signif.", signif_label))) +
    geom_vline(xintercept = c(-logFC, logFC), color = "lightgray", linetype = "dashed") +
    geom_hline(yintercept = -log10(pval), color = "lightgray", linetype = "dashed")
  
  # Add labels if enabled
  if (label) {
    g <- g + ggrepel::geom_text_repel(
      data = topnames,
      aes_string(label = label_column),
      max.overlaps = max.overlaps,
      color = "blue",
      size = label_size,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"),
      show.legend = NA
    )
  }
  
  return(g)
}
