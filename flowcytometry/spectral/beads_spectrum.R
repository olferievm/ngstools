#' Calculate spectral profiles of positive and negative fluorescent beads
#'
#' Identifies positive and negative bead populations using the complete
#' spectral detector matrix and calculates the median fluorescence profile
#' for each population across all supplied detector channels.
#'
#' An optional FSC/SSC or other gate can be applied before population
#' identification. The detector matrix can be transformed with an arcsinh
#' transformation before PCA and clustering. The spectral profiles returned
#' by the function are calculated from the original, untransformed
#' fluorescence values.
#'
#' The function performs PCA on the transformed spectral matrix and uses a
#' two-component Gaussian mixture model to identify the two bead populations.
#' The cluster with the greater overall fluorescence signal is labeled
#' \code{"Pos"}, and the cluster with the lower signal is labeled
#' \code{"Neg"}.
#'
#' @param ff A \code{flowFrame} containing a single fluorescent bead control.
#' @param gate Optional gate to subset \code{ff} before clustering. Typically
#'   an FSC/SSC gate used to remove electronic noise and debris. Must be
#'   compatible with \code{\link[flowCore]{Subset}}.
#' @param channels Character vector specifying the detector channels to use
#'   for spectral analysis. These should correspond to the detector columns
#'   of the spectral matrix.
#' @param transform Logical indicating whether the fluorescence values should
#'   be arcsinh transformed before PCA and clustering. Default is \code{TRUE}.
#' @param cofactor Numeric cofactor used for the arcsinh transformation.
#'   The transformation is \code{asinh(x / cofactor)}. Default is \code{150}.
#' @param npcs Integer specifying the number of principal components used for
#'   clustering. Default is \code{5}.
#'
#' @return A data frame with one row for each detector and two columns,
#'   \code{Neg} and \code{Pos}, containing the median fluorescence intensity 
#'   of the negative and positive bead populations, respectively. The
#'   fluorescence values are calculated from the original, untransformed
#'   expression matrix.
#'
#' @details
#' The function first applies the optional gate, extracts the specified
#' detector channels, and optionally applies an arcsinh transformation.
#' PCA is then performed using the complete detector matrix, and exactly two
#' populations are identified using \code{\link[mclust]{Mclust}}.
#'
#' After clustering, the median fluorescence intensity is calculated
#' independently for each detector using the original fluorescence values.
#' The two clusters are labeled according to their overall fluorescence
#' signal, with the lower-signal cluster designated as \code{Neg} and the
#' higher-signal cluster as \code{Pos}.
#'
#' @examples
#' \dontrun{
#' detectors <- colnames(sp_mat)
#'
#' profile <- beads_spectrum(
#'   ff = fs[["CD56_BUV661.fcs"]],
#'   gate = bead_gate,
#'   channels = detectors
#' )
#'
#' head(profile)
#' }
#'
#' @importFrom flowCore Subset exprs
#' @importFrom matrixStats colMedians
#' @importFrom mclust Mclust
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr %>%
#'
#' @export


beads_spectrum <- function(
    ff,
    gate = NULL,
    channels,
    transform = TRUE,
    cofactor = 150,
    npcs = 5
) {
  
  # FSC/SSC bead gate
  if (!is.null(gate)) {
    ff <- Subset(ff, gate)
  }
  
  # Original fluorescence matrix
  x_raw <- exprs(ff)[, channels, drop = FALSE]
  
  # Matrix used for clustering
  x_cluster <- x_raw
  
  if (transform) {
    x_cluster <- asinh(x_cluster / cofactor)
  }
  
  # PCA using all spectral channels
  xpca <- prcomp(
    x_cluster,
    center = TRUE,
    scale. = TRUE
  )$x
  
  # Two bead populations
  cluster <- mclust::Mclust(
    xpca[, seq_len(min(npcs, ncol(xpca)))],
    G = 2
  )$classification
  
  # Median spectral profile in ORIGINAL fluorescence space
  cluster_medians <- sapply(
    1:2,
    function(k) {
      matrixStats::colMedians(
        x_raw[cluster == k, , drop = FALSE],
        na.rm = TRUE
      )
    }
  )
  
  # Determine which cluster is negative/positive
  cluster_signal <- colSums(cluster_medians)
  
  clust_labels <- c("Neg", "Pos")[order(cluster_signal)]
  
  colnames(cluster_medians) <- clust_labels
  
  # Return spectral profiles
  cluster_medians %>%
    as.data.frame() %>%
    tibble::rownames_to_column(
      var = "Detector"
    )
}