plot_flowset_pages <- function(
    fs,
    plot_fun,
    nrow = 2,
    ncol = 2,
    file,
    width = 11,
    height = 8.5
) {
  
  n_per_page <- nrow * ncol
  
  cairo_pdf(
    file,
    width = width,
    height = height
  )
  
  on.exit(dev.off(), add = TRUE)
  
  for (i in seq(1, length(fs), by = n_per_page)) {
    
    idx <- i:min(i + n_per_page - 1, length(fs))
    fs_page <- fs[idx]
    
    # n <- 2
    n <- length(fs_page)

    # Adjust layout for the last page
    if (n < n_per_page) {
      page_nrow <- ceiling(n / ncol)
      page_ncol <- min(n, ncol)
    } else {
      page_nrow <- nrow
      page_ncol <- ncol
    }
    
    g <- plot_fun(
      fs_page,
      nrow = page_nrow,
      ncol = page_ncol
    )
    
    print(g)
  }
}