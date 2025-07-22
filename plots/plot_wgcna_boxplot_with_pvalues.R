
#' Plot WGCNA Module Eigengene Boxplot with Statistical Annotation
#'
#' This function creates a boxplot comparing WGCNA module eigengene values
#' across timepoints and study groups, with p-values overlaid from a provided
#' model result data frame.
#'
#' @param eigengene_col Character. Name of the column in `tmp` representing the WGCNA eigengene (e.g., "MEpurple").
#' @param term_col Character. Name of the binary annotation column (e.g., "AVE.ACUTEBLOODLOSS_ANEMIA").
#' @param term_label Character. Display label for the annotation used in the legend.
#' @parma ylab Character. Add functional cluster name, aka B cell, to the y axis label.
#' @param xlab Character. x-axis label. Will be shown unmodified.
#' @param color_palette = c("No" = "#377EB8", "Yes" = "#FF7F00") Named vector of colors for ploted groups.
#' @param coef_outlier - NULL (ignored) or numeric value (1.5) to remove outliers among y parameter.
#' @param split_by - NULL (ignored) or factor to split plots by groups.
#' @param top_mixed_results Data frame. Statistical model results with p-values and contrasts (e.g., `top.mes.lmer`).
#' @param use_adjust_pvalue = TRUE - Default is TRUE. Rarely we like to depict real p - values
#' @param tmp Data frame. Contains sample-level data with grouping variables and eigengene values.
#' @param image_path Character. Path where to save the output image.
#' 
#' @return Saves a ggplot boxplot with significance brackets as a PNG file.

plot_wgcna_boxplot_with_pvalues <- function(eigengene_col,
                                            term_col,
                                            term_label,
                                            ylab = "",
                                            xlab = "",
                                            color_palette = c("No" = "#377EB8", "Yes" = "#FF7F00"),
                                            coef_outlier = NULL,
                                            split_by = NULL,
                                            top_mixed_results,
                                            use_adjust_pvalue = TRUE,
                                            tmp,
                                            image_path) {
  # ---- Labels and file name ----
  y_label <- sprintf('ME %s (%s)', ylab, gsub("ME","",eigengene_col))
  file_name <- sprintf('wgcna_boxplot_facet_group_%s_%s.png',
                       eigengene_col,
                       term_col)
  
  # ---- Extract relevant p-values for this eigengene and annotation term ----
  pvl.df <- top_mixed_results %>%
    dplyr::filter(!!sym("wgcnacolor") == eigengene_col,
                  !!sym("term") == term_col) %>%
    tidyr::separate(
      col = 'contrast',
      into = c('group1', 'group2'),
      sep = " - ",
      remove = TRUE
    )
  
  # ---- Clean data ---------------------------------------------------------
  tmp <- tmp %>%
    # ❶ drop rows with any NA in key columns
    dplyr::filter(!is.na(!!sym(term_col)), !is.na(!!sym(eigengene_col))) %>%
    droplevels()
  
  # ---- Optional outlier removal ------------------------------------------
  if (!is.null(coef_outlier)) {
    tmp <- tmp %>%
      dplyr::mutate(.is_outlier = rstatix::is_outlier(!!sym(eigengene_col), coef = coef_outlier)) %>%
      dplyr::filter(.is_outlier == FALSE) %>%
      dplyr::select(-.is_outlier)
  }
  
  # ---- y‑axis position for p‑value brackets ------------------------------
  if (is.null(split_by)) {
    ## ❷ no extra facet – use Event only
    y.df <- tmp %>%
      dplyr::group_by(Event) %>%
      dplyr::summarise(y.position = max(!!sym(eigengene_col), na.rm = TRUE) * 1.20,
                       .groups   = "drop")
    join_vars <- "Event"
  } else {
    ## ❸ facet by an extra grouping column
    y.df <- tmp %>%
      dplyr::group_by(!!sym(split_by), Event) %>%
      dplyr::summarise(y.position = max(!!sym(eigengene_col), na.rm = TRUE) * 1.20,
                       .groups   = "drop")
    join_vars <- c(split_by, "Event")
  }
  
  # ---- If use unadjuct p vales.
  if(!use_adjust_pvalue | !any(colnames(pvl.df) == 'p.adj')){
    pvl.df$p.adj <- pvl.df$pval
  }
  
  # ---- Merge p‑values with y‑positions and prettify labels ---------------
  p.df <- pvl.df %>%
    dplyr::left_join(y.df, by = join_vars) %>%
    dplyr::mutate(
      p.label = dplyr::case_when(
        p.adj <= 0.0005 ~ "***",
        p.adj <= 0.005  ~ "**",
        p.adj <= 0.05   ~ "*",
        p.adj <= 0.10   ~ ".",
        TRUE            ~ "n/s"
      ),
      Event = factor(Event, levels = c("PreOp", "POD1", "W6")),
      xmin  = as.integer(Event) - 0.20,
      xmax  = as.integer(Event) + 0.20
    )
  
  
  
  # ---- Generate plot ----
  g <- ggplot(tmp, aes(x = Event, y = !!sym(eigengene_col))) +
    stat_boxplot(aes(fill = !!sym(term_col)),
                 geom = 'errorbar',
                 position = position_dodge(width = 0.8)) +
    geom_boxplot(aes(fill = !!sym(term_col)),
                 outlier.shape = NA,
                 position = position_dodge(width = 0.8)) +
    geom_jitter(aes(fill = !!sym(term_col)),
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
    scale_fill_manual(values = color_palette) +
    labs(x = xlab, y = y_label, fill = term_label) +
    theme_cowplot() +
    theme(legend.position = "top")
  
  if (!is.null(split_by)) {
    # Here the simple ~ !!sym(split_by) fails.
    facet_sym <- rlang::sym(split_by)
    g <- g + facet_wrap(rlang::enquo(facet_sym))
  }
  
  g <- g +
    ggprism::add_pvalue(
      p.df,
      xmin = "xmin",
      xmax = "xmax",
      label = "p.label",
      fontface = "plain",
      fontfamily = "serif",
      label.size = 5,
      y.position = "y.position",
      tip.length = 0.03
    )
  
  # ---- Save to file ----
  ggplot2::ggsave(
    filename = here::here(image_path, file_name),
    plot = g,
    width = 10,
    height = 5
  )
  
  invisible(g)  # Optionally return plot invisibly
}

