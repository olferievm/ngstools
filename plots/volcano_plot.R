#' Volcano Plot for Gene Expression Results
#'
#' This function generates a volcano plot from a `TopTags` object (from edgeR)
#' or a data frame containing gene expression results with `logFC`, `FDR`
#' (`adj.P.Value`, `adj.pval`), `PValue` (`P.Value`, `pval`), and a gene label column.
#'
#' @param top A `TopTags` object or a data frame with columns `logFC`, `FDR`
#' (`adj.P.Value`, `adj.pval`), `PValue` (`P.Value`, `pval`) and `label_column`.
#' @param logFoldChange Numeric. Threshold for absolute log2 fold change. Default: 1.
#' @param pval Numeric. P value or FDR cutoff significance (between 0 and 1). Default: 0.05
#' @param plotAdjustedPvalue Logical. Plot adjusted p-value (e.g. FDR) instead of raw p-values. Default: FALSE.
#' @param useAdjustedPvalueCutoff Logical. Use adjusted p-value (e.g. FDR) cutoff to color/label genes. Default: TRUE.
#' @param label Logical. Whether to label significant genes. Default: TRUE.
#' @param main Character. Plot title. Default: "Volcano plot".
#' @param xlim "auto" or numeric vector of length 2 for x-axis limits. Default: "auto".
#' @param ylim "auto" or numeric vector of length 2 for y-axis limits. Default: "auto".
#' @param max.overlaps Integer. Max overlapping labels (ggrepel). Default: 10.
#' @param pval_colname Character. Column name for raw p-values. Default: "PValue".
#' @param adjpval_colname Character. Column name for adjusted p-values. Default: "FDR".
#' @param logFC_colname Character. Column name for logFC. Default: "logFC".
#' @param label_column Character. Column name for text labels. Default: "gene_name".
#' @param color_scale Character vector of length 2. Colors for
#'   c("Not Signif.", significant_label). Default: c("black", "orange").
#' @param point_size Numeric. Point size. Default: 1.75.
#' @param point_alpha Numeric. Point opacity. Default: 0.5.
#' @param label_size Numeric. Text size for gene labels. Default: 3.
#' @param max_points Numeric. Downsample background points to this many (keep all significant). Default: 0 (no downsampling).
#'
#' @return A `ggplot2` object.
#' @import ggplot2
#' @import ggrepel
#' @import dplyr
#' @export
volcano_plot <- function(top,
                         logFoldChange = 1,
                         pval = 0.05,
                         plotAdjustedPvalue = FALSE,
                         useAdjustedPvalueCutoff = TRUE,
                         label = TRUE,
                         main = "Volcano plot",
                         xlim = "auto",
                         ylim = "auto",
                         max.overlaps = 10,
                         pval_colname = "PValue",
                         adjpval_colname = "FDR",
                         logFC_colname = "logFC",
                         label_column = "gene_name",
                         color_scale = c("black", "orange"),
                         point_size = 1.75,
                         point_alpha = 0.5,
                         label_size = 3,
                         max_points = 0) {
  
  # Convert TopTags
  if (inherits(top, "TopTags")) top <- top$table
  if (!is.data.frame(top)) stop("'top' must be a data.frame or a TopTags object.")
  
  # Check required columns
  required_cols <- c(logFC_colname, pval_colname, adjpval_colname, label_column)
  missing_cols <- setdiff(required_cols, colnames(top))
  if (length(missing_cols) > 0) {
    stop("Missing column(s): ", paste(missing_cols, collapse = ", "),
         "\nCheck input or adjust: pval_colname, adjpval_colname, logFC_colname, label_column")
  }
  
  # Warn if inconsistent cutoff/plot options
  if (plotAdjustedPvalue && !useAdjustedPvalueCutoff) {
    warning("Using adjusted p-values for plotting but unadjusted cutoff for significance. ",
            "Switching to adjusted cutoff.")
    useAdjustedPvalueCutoff <- TRUE
  }
  
  # Select & rename standard columns
  top <- top %>%
    dplyr::select(all_of(c(logFC_colname, pval_colname, adjpval_colname, label_column)))
  if (logFC_colname != "logFC") top <- top %>% rename(logFC = !!sym(logFC_colname))
  if (pval_colname != "PValue") top <- top %>% rename(PValue = !!sym(pval_colname))
  if (adjpval_colname != "FDR") top <- top %>% rename(FDR = !!sym(adjpval_colname))
  
  # Column logic
  pval_cutoff_colname <- ifelse(useAdjustedPvalueCutoff, "FDR", "PValue")
  pval_plot_colname   <- ifelse(plotAdjustedPvalue, "FDR", "PValue")
  
  # Validity checks
  if (!(logFoldChange > 0)) stop("logFoldChange must be > 0.")
  if (pval <= 0 || pval >= 1) stop("pval cutoff must be between 0 and 1.")
  if (!is.character(color_scale) || length(color_scale) != 2) {
    stop("color_scale must be a character vector of length 2.")
  }
  
  # Axis limits
  if (is.character(xlim) && xlim == "auto") {
    xmax <- ceiling(max(abs(top$logFC), na.rm = TRUE)) + 1
    xlim <- c(-xmax, xmax)
  }
  if (is.character(ylim) && ylim == "auto") {
    actualpval <- ifelse(plotAdjustedPvalue, top$FDR, top$PValue)
    ymax <- ceiling(max(-log10(actualpval), na.rm = TRUE)) + 1
    ylim <- c(0, ymax)
  }
  
  # Significance label
  signif_label <- case_when(
    pval_cutoff_colname == "FDR"   ~ sprintf("FDR < %.2g", pval),
    pval_cutoff_colname == "PValue"~ sprintf("pval < %.2g", pval),
    TRUE                           ~ sprintf("p < %.2g", pval)
  )
  
  # Annotate significance
  top <- top %>%
    mutate(Threshold = ifelse(abs(logFC) > logFoldChange &
                                (!!sym(pval_cutoff_colname)) < pval,
                              signif_label, "Not Signif."))
  
  # Subset sig genes for labels
  topnames <- filter(top, Threshold == signif_label)
  
  # Optional downsampling of background
  if (is.numeric(max_points) && max_points > 0 && max_points < nrow(top)) {
    top$pick <- FALSE
    top$pick[sample(1:nrow(top), max_points)] <- TRUE
    top <- filter(top, Threshold == signif_label | pick)
  }
  
  # y-axis label
  y_label <- if (plotAdjustedPvalue) {
    if (pval_plot_colname == "FDR") "FDR" else "adjPvalue"
  } else {
    "pvalue"
  }
  ylab_expr <- bquote(-log[10](.(y_label)))
  
  # Plot
  g <- ggplot(top, aes(x = logFC,
                       y = -log10(!!sym(pval_plot_colname)),
                       color = Threshold)) +
    geom_point(alpha = point_alpha, size = point_size) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom") +
    xlim(xlim) + ylim(ylim) +
    labs(title = main,
         x = "log2 fold change",
         y = ylab_expr) +
    scale_color_manual(values = setNames(color_scale, c("Not Signif.", signif_label))) +
    geom_vline(xintercept = c(-logFoldChange, logFoldChange),
               color = "lightgray", linetype = "dashed") +
    geom_hline(yintercept = -log10(pval),
               color = "lightgray", linetype = "dashed")
  
  # Labels
  if (label) {
    g <- g + ggrepel::geom_text_repel(
      data = topnames,
      aes(x = logFC,
          y = -log10(!!sym(pval_plot_colname)),
          label = !!sym(label_column)),
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
