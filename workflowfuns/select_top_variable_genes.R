#' Select top variable genes from an expression matrix
#'
#' Identify the most variable genes across samples and return them in a format
#' suitable for downstream network analysis (e.g., WGCNA). Variability can be
#' quantified using variance or median absolute deviation (MAD).
#'
#' @param expr_matrix Numeric matrix of gene expression values with
#'   genes in rows and samples in columns.
#' @param n_top Integer. Number of most variable genes to retain.
#'   Default is 1000.
#' @param method Character. Method used to quantify variability.
#'   Options are `"var"` (variance) or `"mad"` (median absolute deviation).
#'   Default is `"mad"`.
#'
#' @return A numeric matrix containing the selected genes, transposed to
#'   **samples × genes** format (required by WGCNA-style workflows).
#'
#' @details
#' The function performs the following steps:
#' 1. Removes genes with zero variance across samples.
#' 2. Computes a variability score for each gene using either variance or MAD.
#' 3. Ranks genes by variability.
#' 4. Returns the top `n_top` genes or `top` quantile genes. E.g. if top = 0.7 return more than 1000 genes only top 1000 will be kept.
#'
#' MAD is often preferred because it is more robust to outliers than variance.
#'
#' @examples
#' expr <- matrix(rnorm(10000), nrow = 1000)
#' rownames(expr) <- paste0("Gene", 1:1000)
#'
#' top_expr <- select_top_variable_genes(expr, n_top = 500, method = "mad")
#' dim(top_expr)
#'
#' @export
select_top_variable_genes <- function(expr_matrix, top = 0.7, n_top = 1000, method = "mad") {
  
  if (!is.matrix(expr_matrix)) {
    stop("Abort: provide an expression matrix with genes in rows and samples in columns.")
  }
  
  if (nrow(expr_matrix) < 50) {
    warning("Matrix has fewer than 50 genes.")
  }
  
  if (!is.numeric(top) || is.null(top) || length(top) != 1) top <- 0
  top <- max(min(top, 0.95), 0)
  
  if (!is.numeric(n_top) || is.null(n_top) || length(n_top) != 1) {
    n_top <- nrow(expr_matrix)
  }
  
  if (n_top < 2 || n_top > nrow(expr_matrix)) {
    n_top <- nrow(expr_matrix)
  }
  
  # remove zero-variance genes
  keep <- apply(expr_matrix, 1, var, na.rm = TRUE) > 0
  expr_matrix <- expr_matrix[keep, , drop = FALSE]
  
  if (method == "var") {
    score <- apply(expr_matrix, 1, var)
  } else if (method == "mad") {
    score <- apply(expr_matrix, 1, mad)
  } else {
    stop("method must be 'var' or 'mad'")
  }
  
  score <- sort(score, decreasing = TRUE)
  
  if (top > 0) {
    thr <- quantile(score, top, na.rm = TRUE)
    n_max <- sum(score > thr)
  } else {
    n_max <- length(score)
  }
  
  limit <- min(n_top, n_max, length(score))
  top_genes <- names(score)[seq_len(limit)]
  
  # WGCNA format: samples x genes
  t(expr_matrix[top_genes, , drop = FALSE])
}