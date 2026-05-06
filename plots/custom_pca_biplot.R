#' Custom PCA biplot with independent aesthetics for individuals, variables, and supplementary quantitative variables
#'
#' Generate a highly customizable PCA biplot from a \code{FactoMineR::PCA} object.
#' This function provides full manual control over aesthetics for individuals
#' (samples), variables (genes), and optional supplementary quantitative variables
#' (e.g., cytokines), using independent color scales via \code{ggnewscale}.
#'
#' Variable coordinates are internally rescaled to match the individual coordinate
#' space for proper biplot representation. Supplementary quantitative variables
#' (if present in \code{pca$quanti.sup}) can be optionally projected and visualized.
#'
#' @param pca A \code{FactoMineR::PCA} object.
#' @param ind_groups Optional vector defining group assignment for individuals.
#' Length must match number of samples and be named or in the same order.
#' @param ind_color Named vector of colors for individual groups.
#' @param ind_size Numeric. Point size for individuals.
#' @param ind_alpha Numeric. Transparency for individuals.
#' @param ind_shape Optional vector assigning shape per individual.
#' @param shape_format Named vector mapping group levels to shapes.
#' @param var_color Named vector of colors for variable groups.
#' @param var_size Numeric. Line width for variable arrows.
#' @param var_line_shape Line type for variable arrows (e.g., "solid", "dashed").
#' @param var_alpha Numeric. Transparency for variable arrows.
#' @param var_groups Named vector mapping variable (gene) names to groups.
#' @param repel_ind Logical. Use \code{ggrepel} for individual labels.
#' @param repel_ind_size Numeric. Label size for individuals.
#' @param repel_ind_max_overlaps Maximum overlaps for individual labels.
#' @param repel_var Logical. Use \code{ggrepel} for variable labels.
#' @param repel_var_size Numeric. Label size for variables.
#' @param repel_var_max_overlaps Maximum overlaps for variable labels.
#' @param show_quanti_sup Character vector of supplementary quantitative variable
#' names to display. Must exist in \code{pca$quanti.sup$coord}.
#' @param quanti_sup_colors Named vector of colors for supplementary variables.
#' @param quanti_sup_linewidth Numeric. Line width for supplementary arrows.
#' @param quanti_sup_linetype Line type for supplementary arrows.
#' @param quanti_sup_alpha Numeric. Transparency for supplementary arrows.
#' @param repel_quanti Logical. Use \code{ggrepel} for supplementary labels.
#' @param repel_quanti_max_overlaps Maximum overlaps for supplementary labels.
#' @param repel_quanti_size Numeric. Label size for supplementary variables.
#'
#' @details
#' The function rescales variable loadings to match the scale of individual scores
#' using a ratio derived from maximum absolute values along each principal component.
#' This ensures correct geometric representation in the biplot.
#'
#' Color scales for individuals, variables, and supplementary variables are handled
#' independently using \code{ggnewscale::new_scale_color()}.
#'
#' Arrow lengths are determined by PCA coordinates; arrowhead size is fixed and
#' purely aesthetic.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \dontrun{
#' g <- custom_pca_biplot(
#'   pca,
#'   ind_groups = tmp$group,
#'   ind_color = c(FH = "#4DAF4A", FS = "#FF7F00",
#'                 MH = "#377EB8", MS = "#E41A1C"),
#'   ind_shape = tmp$group,
#'   shape_format = c(FH = 16, FS = 16, MH = 17, MS = 17),
#'   var_color = c(IFNa="#b80422", IFNg="#172767"),
#'   var_groups = var_groups,
#'   repel_ind = FALSE,
#'   repel_var = TRUE,
#'   repel_var_size = 3,
#'   repel_var_max_overlaps = Inf,
#'   show_quanti_sup = c("MSD_CTKNS_IFNγ", "IFNA"),
#'   quanti_sup_colors = c("MSD_CTKNS_IFNγ"="#c3d878",
#'                         "IFNA"="#58a787"),
#'   repel_quanti = TRUE,
#'   repel_quanti_max_overlaps = Inf,
#'   repel_quanti_size = 4
#' )
#' }
#'
#' @import ggplot2
#' @import ggrepel
#' @import ggnewscale
#' @export



custom_pca_biplot <- function(
    pca,
    axes = c(1,2),
    ind_groups = NULL,
    ind_color = NULL,
    ind_size = 2,
    ind_alpha = 1,
    ind_shape = NULL,
    shape_format = NULL,
    var_color = NULL,
    var_size = 0.5,
    var_line_shape = 'solid',
    var_alpha = 0.5,
    var_groups = NULL,
    var_top = NULL,
    repel_ind = FALSE,
    repel_ind_size = 2,
    repel_ind_max_overlaps = 50,
    repel_var = TRUE,
    repel_var_size = 3,
    repel_var_max_overlaps = 50,
    show_quanti_sup = NULL,
    quanti_sup_colors = NULL,
    quanti_sup_linewidth = 2,
    quanti_sup_linetype = 'solid',
    quanti_sup_alpha = 0.7,
    repel_quanti = TRUE,
    repel_quanti_max_overlaps = Inf,
    repel_quanti_size = 3
) {
  library(ggplot2)
  library(ggrepel)
  library(ggnewscale)
  library(cowplot)
  
  
  # Contribution of axes
   rate_of_variance <- signif(pca$eig[axes,"percentage of variance"],2)
   xaxis <- sprintf("Dim.%s",axes[1])
   yaxis <- sprintf("Dim.%s",axes[2])
   xlab <- sprintf("Dim%02s %s%%",axes[1],rate_of_variance[1])
   ylab <- sprintf("Dim%02s %s%%",axes[2],rate_of_variance[2])
  
   # --- Calculate scale factors coordinates ---
   scale_factor_ind <- apply(pca$ind$coord,2,function(x){max(abs(x))})
   scale_factor_var <- apply(pca$var$coord,2,function(x){max(abs(x))})
   scale_ratio <- 0.8 * (scale_factor_ind / scale_factor_var)
   
  # --- Extract coordinates individuals ---
  ind_df <- as.data.frame(pca$ind$coord)
  ind_df$sample <- rownames(ind_df)
  
  # --- Extract coordinates of vaiables and scale them ---
  # Multiply var$coord by a multiplier that is typically the ratio
  # of standard deviations or the maximum value of observations
  #  divided by the maximum value of variables
  var_df <- as.data.frame(pca$var$coord * scale_ratio)
  var_df$gene <- rownames(var_df)
  
  # --- Individual grouping ---
  if (!is.null(ind_groups)) {
    ind_df$group <- ind_groups
  } else {
    ind_df$group <- "all"
  }
  
  # --- Individual shape ---
  if (!is.null(ind_shape)) {
    ind_df$shape <- ind_shape
  } else {
    ind_df$shape <- "all"
  }
  
  # --- Variable grouping ---
  if (!is.null(var_groups)) {
    var_df$group <- var_groups[var_df$gene]
  }else{
    var_df$group <- "other"
  }
  
  # Select only genes with top contrib (this should be redone)
  if(!is.null(var_top)){
    top_axes_1 <- names(sort(pca$var$contrib[,axes[1]], decreasing = TRUE)[1:var_top])
    top_axes_2 <- names(sort(pca$var$contrib[,axes[2]], decreasing = TRUE)[1:var_top])
    sel_vars <- unique(c(top_axes_1,top_axes_2))
    var_df <- var_df[sel_vars,]
  }
  
  if(!is.null(show_quanti_sup) || !is.null(pca$quanti.sup)){
    quanti_sup <- pca$quanti.sup$coord[show_quanti_sup,]
    quanti_sup <- quanti_sup * scale_ratio
    quanti_sup <- as.data.frame(quanti_sup)
    quanti_sup$quanti <- rownames(quanti_sup)
  }else{
    quanti_sup <- NULL
  }
  
  # --- Base ---
  g <- ggplot()
  
  # =========================
  # Individuals (first scale)
  # =========================
  g <- g + geom_point(
    data = ind_df,
    aes(x = !!sym(xaxis), y = !!sym(yaxis),
        color = group,
        shape = shape),
    alpha = ind_alpha,
    size = ind_size
  )

  if (repel_ind) {
    g <- g + geom_text_repel(
      data = ind_df,
      aes(x = !!sym(xaxis), y = !!sym(yaxis), label = sample, color = group),
      max.overlaps = repel_ind_max_overlaps,
      size = repel_ind_size
    )
  }
  
  
  if (!is.null(ind_color)) {
    g <- g + scale_color_manual(values = ind_color)
  }
  
  if (!is.null(ind_shape)) {
    g <- g + scale_shape_manual(values = shape_format)
  }
  
  # reset color scale
  g <- g + ggnewscale::new_scale_color()
  
  # =========================
  # Variables (second scale)
  # =========================
  g <- g + geom_segment(
    data = var_df,
    aes(x = 0, y = 0,
        xend = !!sym(xaxis), yend = !!sym(yaxis),
        color = group),
    alpha = var_alpha,
    linewidth = var_size,
    linetype = var_line_shape,
    arrow = grid::arrow(length = grid::unit(0.15, "cm"))
  )
  
  if (repel_var) {
    g <- g + geom_text_repel(
      data = var_df,
      aes(x = !!sym(xaxis), y = !!sym(yaxis), label = gene, color = group),
      size = repel_var_size,
      max.overlaps = repel_var_max_overlaps
    )
  }
  
  if (!is.null(var_color)) {
    g <- g + scale_color_manual(values = var_color)
  }
  
  if (!is.null(quanti_sup)){

  # reset color scale
  g <- g + ggnewscale::new_scale_color()
  g <- g + geom_segment(
    data = quanti_sup,
    aes(x = 0, y = 0,
        xend = !!sym(xaxis), yend = !!sym(yaxis),
        color = quanti),
    linewidth = quanti_sup_linewidth,
    linetype = quanti_sup_linetype,
    arrow = grid::arrow(length = grid::unit(0.15, "cm"))
  )
  if(!is.null(quanti_sup_colors)){
    g <- g + scale_color_manual(values = quanti_sup_colors)
  }
  if(repel_quanti){
    g <- g + geom_text_repel(
      data = quanti_sup,
      aes(x = !!sym(xaxis), y = !!sym(yaxis), label = quanti, color = quanti),
      size = repel_quanti_size,
      max.overlaps = repel_quanti_max_overlaps
    )}
  }
  
  g <- g + theme_cowplot() + xlab(xlab) + ylab(ylab)
  return(g)
}