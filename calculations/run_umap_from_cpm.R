#' Run PCA → UMAP on a CPM Matrix
#'
#' This function takes a CPM expression matrix, performs log-transformation,
#' principal component analysis (PCA), and then runs UMAP on the first
#' `n_pcs` principal components. The function automatically ensures that samples
#' are rows and genes are columns.
#'
#' @param cpm_mat A numeric matrix or data frame of CPM values.  
#'   If genes > samples, the matrix will be transposed so that rows = samples.
#'
#' @param n_pcs Integer. The number of principal components to use for UMAP  
#'   (recommended: 10–15). Default = 15.
#'
#' @param n_neighbors Integer. Number of nearest neighbors for UMAP.  
#'   Controls local vs global structure. Default = 30.
#'
#' @param min_dist Numeric. Controls compactness of UMAP clusters  
#'   (lower = tighter clusters). Default = 0.3.
#'
#' @param seed Integer. Random seed for reproducible PCA/UMAP results.  
#'   Default = 123.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{pca}{A `prcomp` object containing PCA results.}
#'     \item{pcs}{A matrix of the first `n_pcs` principal components.}
#'     \item{umap}{A data frame containing UMAP1, UMAP2, and sample names.}
#'   }
#'
#' @details
#' The function applies a log2(CPM + 1) transformation before PCA.  
#' UMAP is computed using the `uwot` package with Euclidean distance.
#'
#' @examples
#' \dontrun{
#' cpm <- read.csv("cpm_matrix.csv", row.names = 1)
#' result <- run_umap_from_cpm(cpm, n_pcs = 12)
#'
#' library(ggplot2)
#' ggplot(result$umap, aes(UMAP1, UMAP2)) +
#'   geom_point() +
#'   theme_minimal()
#' }
#'
#' @export
#' 

run_umap_from_cpm <- function(cpm_mat, n_pcs = 15, n_neighbors = 30, min_dist = 0.3,
                              seed = 123) {
  set.seed(seed)

  require(uwot)
  require(scales)
  
  # Ensure samples are rows
  if (nrow(cpm_mat) > ncol(cpm_mat)) {
    # assume genes > samples → transpose
    message("Transposing input so that rows = samples, columns = genes.")
    cpm_mat <- t(cpm_mat)
  }
  
  # Log-transform CPM values
  log_cpm <- log2(cpm_mat + 1)
  
  ## ---- PCA ----
  pca <- prcomp(log_cpm, center = TRUE, scale. = TRUE)
  
  # Select first n PCs (recommended 10–15)
  pcs <- pca$x[, 1:n_pcs, drop = FALSE]
  
  ## ---- UMAP ----
  umap_res <- umap(
    pcs, 
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = "euclidean",
    n_threads = 4,
    ret_model = FALSE
  )
  
  umap_df <- data.frame(
    UMAP1 = umap_res[,1],
    UMAP2 = umap_res[,2],
    sample = rownames(umap_res)
  )
  
  ## ---- Return results ----
  return(list(
    pca = pca,
    pcs = pcs,
    umap = umap_df
  ))
}

