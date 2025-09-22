#' Volcano Plot for Gene Expression Results
#'
#' This function generates a volcano plot from a `TopTags` object (from edgeR)
#' or a data frame containing gene expression results with `logFC`, `FDR` (`adj.P.Value`, `adj.pval`), `PValue` (`P.Value`,`pval`) and
#' a gene label column.
#'
#' @param top A `TopTags` object or a data frame with columns `logFC`, `FDR` (`adj.P.Value`, `adj.pval`, `PValue`, `P.Value`, `pval`) and `label_column`.
#' @param logFC Numeric. Threshold for absolute log2 fold change. Default: 1.
#' @param pval Numeric. FDR cutoff for significance (between 0 and 1). Default: 0.05.
#' @param label Logical. Whether to label significant genes. Default: TRUE.
#' @param main Character. Plot title. Default: "Volcano plot".
#' @param xlim "auto" or numeric vector of length 2 for x-axis limits. Default: "auto".
#' @param ylim "auto" or numeric vector of length 2 for y-axis limits. Default: "auto".
#' @param max.overlaps Integer. Max overlapping labels (ggrepel). Default: 10.
#' @param pval_colname Character. Column name to use for significance. Default: "FDR".
#' @param logFC_colname Character. Column name to use for logFC. Default: "logFC".
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
                         pval_colname = 'FDR',
                         logFC_colname = 'logFC',
                         label_column = "gene_name",
                         color_scale = c("black", "orange"),
                         label_size = 3,
                         max_points = FALSE) {
  
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
  required_cols <- c(logFC_colname, pval_colname, label_column)
  missing_cols <- setdiff(required_cols, colnames(top))
  
  if (length(missing_cols) > 0)
    stop("Missing column(s): ", paste(missing_cols, collapse = ", "), "\n - Check the input or adjust: pval_column, logFC_column, label_column")
  
  top <- top %>% dplyr::select(all_of(c(logFC_colname, pval_colname, label_column)))
  
  if(logFC_colname != 'logFC'){
    top <- top %>% dplyr::rename(FDR = !!sym(logFC_colname))
  }
  if(pval_colname != 'FDR'){
    top <- top %>% dplyr::rename(FDR = !!sym(pval_colname))
  }
  
  
  # Check input validity
  if (!(logFC > 0)) stop("logFC must be > 0.")
  if (pval <= 0 || pval >= 1) stop("pval cut-off must be between 0 and 1.")
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
  signif_label = case_when(
    pval_colname == 'FDR' ~ sprintf("FDR < %.2g", pval),
    pval_colname == 'PValue' ~ sprintf("pval < %.2g", pval),
    pval_colname == 'P.Value' ~ sprintf("pval < %.2g", pval),
    pval_colname == 'pval' ~ sprintf("pval < %.2g", pval),
    pval_colname == 'adj.p.value' ~ sprintf("adj.pval < %.2g", pval),
    pval_colname == 'adj.pval' ~ sprintf("adj.pval < %.2g", pval),
    pval_colname == 'apval' ~ sprintf("adj.pval < %.2g", pval),
  )
  
  # Annotate data with significance
  top$Threshold <- ifelse(abs(top$logFC) > logFC & top$FDR < pval,
                          signif_label, "Not Signif.")
  
  # Subset significant points for labeling
  topnames <- dplyr::filter(top, Threshold == signif_label)
  
  # Plot size reduction. Reduced unlabeled points.
  if(max_points > 0){
    if(max_points < nrow(top)){
       top$pick <- FALSE
       top$pick[sample(1:nrow(top), max_points)] <- TRUE
       top <- top %>% dplyr::filter(Threshold == signif_label | pick)
    }
  }
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
      aes(label = !!sym(label_column)),
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
