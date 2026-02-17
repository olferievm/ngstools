ggplot.2groups.comparisons <- function(
    x,
    name,                  # numeric outcome column
    categories,            # x-axis factor column
    groups,                # 2-level grouping column
    test = c("wilcox", "t_test"),
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
    coord_flip = TRUE
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
  
  if (!keep.nonsignif)
    stat.test <- dplyr::filter(stat.test, p.adj.signif != "ns")
  
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