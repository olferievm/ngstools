#' Boxplot comparison of two groups across multiple categories
#'
#' Generate a ggplot2 boxplot with jittered observations and statistical
#' comparisons between two groups performed independently within each category.
#' The function performs either a Wilcoxon rank-sum test or Student's t-test
#' for each category, adjusts p-values using the Benjamini–Hochberg procedure,
#' and displays significance brackets on the plot.
#'
#' Optionally, the function can return only the calculated statistics table
#' without generating a plot, or use externally supplied p-values (e.g.
#' moderated statistics from limma).
#'
#' @param x A data.frame containing the variables to be plotted.
#'
#' @param name Character string specifying the numeric outcome variable.
#'
#' @param categories Character string specifying the categorical variable
#' defining the x-axis categories. Each category will be tested independently.
#'
#' @param groups Character string specifying the grouping variable. Must have
#' exactly two levels.
#'
#' @param test Statistical test used for comparisons. Either `"wilcox"`
#' (Wilcoxon rank-sum test) or `"t_test"` (Student's t-test).
#'
#' @param only_signifs Logical. If `TRUE`, the function returns the table of
#' calculated statistics and does not generate a plot.
#'
#' @param manual_signifs Optional data.frame containing externally calculated
#' significance values. Must contain the same column name as `categories`
#' and the columns `p`, `p.adj`, and `p.adj.signif`. When provided, these
#' values replace internally computed statistics.
#'
#' @param xlab Optional x-axis label.
#'
#' @param ylab Optional y-axis label.
#'
#' @param legend.color Legend title for point colors.
#'
#' @param legend.fill Legend title for boxplot fill colors.
#'
#' @param point.colors Vector of colors used for points corresponding to
#' the levels of `groups`.
#'
#' @param point.opacity Numeric transparency value for jittered points.
#'
#' @param point.size Numeric size of jittered points.
#'
#' @param jitter.width Width of horizontal jitter applied to points.
#'
#' @param dodge.width Dodge width used to separate the two groups.
#'
#' @param fill.colors Vector of fill colors for boxplots corresponding to
#' the levels of `groups`.
#'
#' @param staplewidth Width of boxplot staples.
#'
#' @param keep.nonsignif Logical. If `FALSE`, comparisons with non-significant
#' adjusted p-values (`p.adj.signif == "ns"`) are removed from the plot.
#'
#' @param bracket.size Line width of significance brackets.
#'
#' @param tip.length Length of bracket tips.
#'
#' @param coord_flip Logical. If `TRUE`, the coordinate system is flipped
#' for horizontal boxplots.
#'
#' @param verbose Logical. If `TRUE`, status messages are printed during
#' execution.
#'
#' @return
#' If `only_signifs = TRUE`, returns a data.frame containing statistical test
#' results with adjusted p-values and plotting coordinates. Otherwise returns
#' a `ggplot` object.
#'
#' @details
#' For each level of `categories`, the function compares the two levels of
#' `groups` using the specified statistical test. P-values are adjusted
#' using the Benjamini–Hochberg method to control the false discovery rate.
#'
#' Significance brackets are added using `ggprism::add_pvalue`, with
#' positions automatically calculated based on category index and maximum
#' observed value.
#'
#' @examples
#' 
#' \dontrun{ 
#' set.seed(1)
#' 
#' df <- data.frame(
#'   category = rep(c("A","B","C"), each = 40),
#'   group = rep(rep(c("ctrl","case"), each = 20), 3),
#'   value = rnorm(120) + rep(c(0,0.5,1), each=40)
#' )
#' 
#' ggplot.2groups.comparisons(
#'   x = df,
#'   name = "value",
#'   categories = "category",
#'   groups = "group"
#' )
#' }
#' 
#' \dontrun{
#' ggplot.2groups.comparisons(
#'   x = df,
#'   name = "expression",
#'   categories = "celltype",
#'   groups = "condition",
#'   test = "wilcox",
#'   fill.colors = c("#1B9E77", "#D95F02"),
#'   point.colors = c("#1B9E77", "#D95F02"),
#'   xlab = "Cell type",
#'   ylab = "Expression"
#' )
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import rstatix
#' @import ggprism
#' @import tibble
#'
#' @export




ggplot.2groups.comparisons <- function(
    x,
    name,                  # numeric outcome column
    categories,            # x-axis factor column
    groups,                # 2-level grouping column
    test = c("wilcox", "t_test"), # which test to use
    only_signifs = FALSE, # if TRUE do not produce plot but return significance table
    manual_signifs = NULL, # a user provided data.frame of match significance (e.g. limma moderated test) with the same column name as categories and columns: p, p.adj, and p.adj.signif
    xlab = NULL,
    ylab = NULL,
    legend.color = NULL,
    legend.fill  = NULL,
    point.colors = NULL,
    point.opacity = 0.5,
    point.size = 1,
    jitter.width = 0.25,
    dodge.width = 0.8,
    fill.colors = NULL,
    staplewidth = 0.25,
    keep.nonsignif = FALSE,
    bracket.size = 0.3,
    tip.length = 0.01,
    coord_flip = TRUE,
    verbose = FALSE
){
  
  test <- match.arg(test)
  
  ## ---- enforce factors ----
  x[[categories]] <- factor(x[[categories]])
  x[[groups]]     <- factor(x[[groups]])
  
  if (nlevels(x[[groups]]) != 2)
    stop("groups must have exactly two levels")
  
  ## ---- statistical test per category ----
  stat.test <- dplyr::group_by(x, .data[[categories]])
  
  stat.test <- switch(
    test,
    wilcox = rstatix::wilcox_test(stat.test, as.formula(paste(name, "~", groups))),
    t_test = rstatix::t_test(stat.test,     as.formula(paste(name, "~", groups)))
  )
  
  stat.test <- stat.test |>
    rstatix::adjust_pvalue(method = "BH") |>
    rstatix::add_significance("p.adj")
  
  ## ---- compute numeric x positions for brackets ----
  lvl_map <- tibble::tibble(
    !!categories := levels(x[[categories]]),
    xnum = seq_along(levels(x[[categories]]))
  )
  
  stat.test <- dplyr::left_join(stat.test, lvl_map, by = categories) |>
    dplyr::mutate(
      xmin = xnum - dodge.width/2,
      xmax = xnum + dodge.width/2
    )
  
  ## ---- y positions (top of each panel) ----
  y.df <- x |>
    dplyr::group_by(.data[[categories]]) |>
    dplyr::summarise(
      y.position = max(.data[[name]], na.rm = TRUE) * 1.08,
      .groups = "drop"
    )
  
  stat.test <- dplyr::left_join(stat.test, y.df, by = categories)
  if(verbose){cat('statistics was calculated\n')}
  # Return only calculated significance table.
  if(only_signifs){
    if(verbose){cat('return only statistics\n')}
    return(stat.test)}
  
  #
  if(!is.null(manual_signifs)){
    if(length(intersect(colnames(manual_signifs),c(categories,"p", "p.adj", "p.adj.signif")))==4){
      stat.test <- stat.test %>%
          dplyr::select(-p, -p.adj, -p.adj.signif) %>%
          dplyr::left_join(., limma.pvals, by = categories, unmatched = "drop", keep = FALSE)
      if(verbose){cat('substitute statistics with manual entry\n')}
    }
  }
  
  if (!keep.nonsignif){
    stat.test <- dplyr::filter(stat.test, p.adj.signif != "ns")
    if(verbose){cat('removed non-significant values\n')}
  }
  
  ## ---- plot ----
  g <- ggplot2::ggplot(x, ggplot2::aes(.data[[categories]], .data[[name]])) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = .data[[groups]]),
      position = ggplot2::position_dodge(width = dodge.width),
      outlier.shape = NA,
      staplewidth = staplewidth
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = .data[[groups]]),
      alpha = point.opacity,
      size = point.size,
      position = ggplot2::position_jitterdodge(
        jitter.width = jitter.width,
        dodge.width  = dodge.width
      )
    ) +
    ggplot2::scale_fill_manual(values = fill.colors, name = legend.fill) +
    ggplot2::scale_color_manual(values = point.colors, name = legend.color)
  
  if (coord_flip)
    g <- g + ggplot2::coord_flip()
  
  g <- g +
    ggprism::add_pvalue(
      stat.test,
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      label = "p.adj.signif",
      bracket.size = bracket.size,
      tip.length = tip.length
    )
  
  if(!is.null(xlab)){
    g <- g + xlab(xlab)
  }
  
  if(!is.null(ylab)){
    g <- g + ylab(ylab)
  }
  
  g
}

