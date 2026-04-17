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
#' 4. Returns the top `n_top` genes.
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
select_top_variable_genes <- function(expr_matrix, n_top = 1000, method = "mad") {
  
  # remove zero-variance genes
  expr_matrix <- expr_matrix[apply(expr_matrix, 1, var) > 0, ]
  
  if (method == "var") {
    score <- apply(expr_matrix, 1, var)
  } else if (method == "mad") {
    score <- apply(expr_matrix, 1, mad)
  }
  
  score <- sort(score, decreasing = TRUE)
  top_genes <- names(score)[seq_len(min(n_top, length(score)))]
  
  # WGCNA format: samples x genes
  t(expr_matrix[top_genes, , drop = FALSE])
}