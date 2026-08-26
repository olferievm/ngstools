plotHeatMatrix <- function(
    x,
    p = NULL,
    col_palette = NULL,
    scale.type = c("auto", "diverging", "positive", "negative"),
    text.size = 3,
    axis.text.size = 10,
    text.color = "green",
    label.precision = 3,
    label.threshold = 0.05,
    reorder.rows = TRUE,
    reorder.cols = TRUE,
    xmin = NULL,
    xmid = NULL,
    xmax = NULL,
    main = "",
    xlab = "",
    ylab = ""
) {
  
  # ============================================================
  # Input checks
  # ============================================================
  
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("x must be a numeric matrix.")
  }
  
  if (is.null(rownames(x))) {
    rownames(x) <- seq_len(nrow(x))
  }
  
  if (is.null(colnames(x))) {
    colnames(x) <- seq_len(ncol(x))
  }
  
  if (!is.null(p)) {
    
    if (!is.matrix(p) || !is.numeric(p)) {
      stop("p must be a numeric matrix or NULL.")
    }
    
    if (!all(dim(x) == dim(p))) {
      stop("x and p must have the same dimensions.")
    }
  }
  
  scale.type <- match.arg(scale.type)
  
  
  # ============================================================
  # Default palette
  # ============================================================
  
  if (is.null(col_palette)) {
    
    col_palette <- c(
      "#0803A1",
      "#1514A6",
      "#225EA8",
      "#74A9CF",
      "#BDC9E1",
      "#F1EEF6",
      "#FDCC8A",
      "#FC8D59",
      "#D7301F",
      "#C11C38",
      "#B01919"
    )
  }
  
  if (length(col_palette) < 3) {
    stop("col_palette must contain at least 3 colors.")
  }
  
  
  # ============================================================
  # Determine observed data range
  # ============================================================
  
  x_values <- as.numeric(x)
  x_values <- x_values[is.finite(x_values)]
  
  if (length(x_values) == 0) {
    stop("x contains no finite values.")
  }
  
  data_min <- min(x_values)
  data_max <- max(x_values)
  
  
  # ============================================================
  # Automatic scale detection
  # ============================================================
  
  if (scale.type == "auto") {
    
    if (data_min < 0 && data_max > 0) {
      
      scale.type <- "diverging"
      
    } else if (data_min >= 0) {
      
      scale.type <- "positive"
      
    } else {
      
      scale.type <- "negative"
    }
  }
  
  
  # ============================================================
  # Determine palette and limits
  # ============================================================
  
  palette_mid <- ceiling(length(col_palette) / 2)
  
  if (scale.type == "diverging") {
    
    # Full palette:
    #
    # negative → neutral → positive
    
    palette_use <- col_palette
    
    if (is.null(xmin)) {
      xmin <- data_min
    }
    
    if (is.null(xmax)) {
      xmax <- data_max
    }
    
    if (is.null(xmid)) {
      xmid <- 0
    }
    
    if (xmin >= xmid || xmax <= xmid) {
      stop(
        "For scale.type = 'diverging', ",
        "xmin < xmid < xmax is required."
      )
    }
    
    
  } else if (scale.type == "positive") {
    
    # Neutral → positive half
    
    palette_use <- col_palette[palette_mid:length(col_palette)]
    
    if (is.null(xmin)) {
      xmin <- 0
    }
    
    if (is.null(xmax)) {
      xmax <- data_max
    }
    
    if (xmin >= xmax) {
      stop("xmin must be smaller than xmax.")
    }
    
    
  } else if (scale.type == "negative") {
    
    # Negative half → neutral
    
    palette_use <- col_palette[1:palette_mid]
    
    if (is.null(xmin)) {
      xmin <- data_min
    }
    
    if (is.null(xmax)) {
      xmax <- 0
    }
    
    if (xmin >= xmax) {
      stop("xmin must be smaller than xmax.")
    }
  }
  
  
  # ============================================================
  # Prepare data
  # ============================================================
  
  df <- reshape2::melt(
    x,
    varnames = c("row", "column"),
    value.name = "value"
  )
  
  
  # ============================================================
  # Optional significance labels
  # ============================================================
  
  if (!is.null(p)) {
    
    df_p <- reshape2::melt(
      p,
      varnames = c("row", "column"),
      value.name = "p.value"
    )
    
    df <- dplyr::left_join(
      df,
      df_p,
      by = c("row", "column")
    )
    
    pcut <- 10^(-label.precision)
    
    df <- df %>%
      dplyr::mutate(
        labelm = dplyr::case_when(
          
          is.na(p.value) ~ "",
          
          p.value < label.threshold &
            p.value >= pcut ~
            format(
              round(p.value, label.precision),
              nsmall = label.precision
            ),
          
          p.value < pcut ~
            paste0("< ", pcut),
          
          TRUE ~ ""
        )
      )
    
  } else {
    
    df$labelm <- ""
  }
  
  
  # ============================================================
  # Reorder rows
  # ============================================================
  
  if (reorder.rows && nrow(x) > 1) {
    
    h <- hclust(
      dist(x),
      method = "ward.D2"
    )
    
    df$row <- factor(
      df$row,
      levels = rownames(x)[h$order]
    )
  }
  
  
  # ============================================================
  # Reorder columns
  # ============================================================
  
  if (reorder.cols && ncol(x) > 1) {
    
    h <- hclust(
      dist(t(x)),
      method = "ward.D2"
    )
    
    df$column <- factor(
      df$column,
      levels = colnames(x)[h$order]
    )
  }
  
  
  # ============================================================
  # Construct plot
  # ============================================================
  
  g <- ggplot(
    df,
    aes(
      x = column,
      y = row,
      fill = value
    )
  ) +
    geom_tile()
  
  
  # ============================================================
  # Add significance labels only when p was supplied
  # ============================================================
  
  if (!is.null(p)) {
    
    g <- g +
      geom_text(
        data = dplyr::filter(
          df,
          !is.na(labelm),
          labelm != ""
        ),
        aes(label = labelm),
        size = text.size,
        color = text.color
      )
  }
  
  
  # ============================================================
  # Color scale
  # ============================================================
  
  if (scale.type == "diverging") {
    
    g <- g +
      scale_fill_gradientn(
        colours = palette_use,
        limits = c(xmin, xmax),
        rescaler = function(x, ...) {
          scales::rescale_mid(
            x,
            to = c(0, 1),
            from = c(xmin, xmax),
            mid = xmid
          )
        },
        oob = scales::squish
      )
    
  } else {
    
    g <- g +
      scale_fill_gradientn(
        colours = palette_use,
        limits = c(xmin, xmax),
        oob = scales::squish
      )
  }
  
  
  # ============================================================
  # Labels and theme
  # ============================================================
  
  g <- g +
    labs(
      title = main,
      x = xlab,
      y = ylab
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(
        size = axis.text.size
      )
    )
  
  return(g)
}