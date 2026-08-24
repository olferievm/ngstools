#' Generate PCA and UMAP coordinates from spectral flow cytometry data
#'
#' Extracts fluorescence or spectral channels from a \code{flowFrame},
#' optionally applies a gate and arcsinh transformation, performs principal
#' component analysis (PCA), and generates a two-dimensional UMAP embedding.
#'
#' The function is intended for exploratory analysis and clustering of
#' spectral flow cytometry data, including pooled or individually stained
#' bead controls. PCA coordinates can be used directly for downstream
#' clustering, while the UMAP coordinates provide a low-dimensional
#' visualization of the spectral data.
#'
#' @param ff A \code{flowFrame} containing flow cytometry events.
#'
#' @param gate Optional gate used to subset \code{ff} before dimensional
#'   reduction. If \code{NULL}, all events in \code{ff} are retained.
#'
#' @param channels Character vector specifying the fluorescence or spectral
#'   channels to use for transformation, PCA, and UMAP. All supplied channel
#'   names must be present in the expression matrix of \code{ff}.
#'
#' @param transform Logical. If \code{TRUE}, applies an arcsinh
#'   transformation to the selected \code{channels} using \code{cofactor}.
#'   The transformed values are also inserted into the returned expression
#'   matrix. Defaults to \code{TRUE}.
#'
#' @param cofactor Positive numeric value specifying the cofactor used for
#'   the arcsinh transformation:
#'   \deqn{\operatorname{asinh}(x / \mathrm{cofactor})}
#'   Defaults to \code{150}.
#'
#' @param npcs Either \code{"all"} or a positive numeric value specifying
#'   the number of principal components to retain. When \code{"all"}, all
#'   available principal components are returned. PCA is always calculated
#'   using all selected channels before optional subsetting of the returned
#'   coordinates.
#'
#' @param n_neighbors Number of nearest neighbors used by
#'   \code{\link[uwot]{umap}}. Defaults to \code{15}.
#'
#' @param umap_scale Logical indicating whether the input spectral matrix
#'   should be scaled internally by \code{\link[uwot]{umap}} before
#'   constructing the embedding. Defaults to \code{TRUE}.
#'
#' @param seed Optional integer used to set the random seed before UMAP
#'   calculation. If \code{NULL}, no seed is set.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{exprs}}{Expression matrix containing all parameters from the
#'   input \code{flowFrame}. When \code{transform = TRUE}, the selected
#'   \code{channels} contain arcsinh-transformed values. All other parameters
#'   remain unchanged.}
#'
#'   \item{\code{xpca}}{Matrix of PCA coordinates for each event. Rows
#'   correspond to events and columns correspond to retained principal
#'   components.}
#'
#'   \item{\code{umap}}{Two-column matrix containing the UMAP coordinates,
#'   named \code{UMAP1} and \code{UMAP2}.}
#' }
#'
#' @details
#' PCA is performed using the selected fluorescence or spectral channels after
#' optional arcsinh transformation. Variables are centered and scaled before
#' PCA using \code{\link[stats]{prcomp}}.
#'
#' The UMAP embedding is calculated independently from the PCA coordinates,
#' using the transformed spectral matrix directly. Consequently, the returned
#' UMAP is an exploratory representation of the full selected spectral space,
#' whereas \code{xpca} can be used as input for downstream methods such as
#' Gaussian mixture modeling with \code{mclust::Mclust}.
#'
#' If a gate is supplied, gating is performed before extraction of the
#' expression matrix and before all subsequent transformations and dimensional
#' reduction.
#'
#' @seealso
#' \code{\link[flowCore]{Subset}},
#' \code{\link[flowCore]{exprs}},
#' \code{\link[stats]{prcomp}},
#' \code{\link[uwot]{umap}},
#' \code{\link[mclust]{Mclust}}
#'
#' @examples
#' \dontrun{
#' result <- beads_umap(
#'   ff = bead_ff,
#'   channels = spectral_channels,
#'   cofactor = 150,
#'   npcs = "all",
#'   n_neighbors = 15,
#'   seed = 123
#' )
#'
#' # PCA coordinates for clustering
#' xpca <- result$xpca
#'
#' # UMAP coordinates for visualization
#' umap <- result$umap
#'
#' # Gaussian mixture clustering
#' fit <- mclust::Mclust(
#'   xpca,
#'   G = 15:25
#' )
#' }
#'
#' @export

beads_umap <- function(
    ff,
    gate = NULL,
    channels,
    transform = TRUE,
    cofactor = 150,
    npcs = "all",
    n_neighbors = 15,
    umap_scale = TRUE,
    seed = NULL
) {
  
  # ------------------------------------------------------------
  # 1. Apply optional gate
  #
  # Typically used to retain the bead population based on
  # FSC/SSC or another pre-defined gate.
  # ------------------------------------------------------------
  if (!is.null(gate)) {
    ff <- flowCore::Subset(ff, gate)
  }
  
  
  # ------------------------------------------------------------
  # 2. Extract the complete expression matrix
  #
  # 'raw_ff' retains all original parameters. Only the specified
  # fluorescence channels are used for transformation, PCA,
  # clustering, and UMAP.
  # ------------------------------------------------------------
  cat(".")
  raw_ff <- flowCore::exprs(ff)

  # Validate requested fluorescence channels
  missing_channels <- setdiff(channels, colnames(raw_ff))

  if (length(missing_channels) > 0) {
    stop(
      "The following channels are not present in the flowFrame: ",
      paste(missing_channels, collapse = ", ")
    )
  }
  
  
  # ------------------------------------------------------------
  # 3. Extract fluorescence matrix for dimensional reduction
  #
  # Rows = events
  # Columns = fluorescence / spectral channels
  # ------------------------------------------------------------

  x_cluster <- raw_ff[, channels, drop = FALSE]

  
  # ------------------------------------------------------------
  # 4. Optional arcsinh transformation
  #
  # This is applied only to fluorescence channels used for PCA
  # and UMAP. The transformed values are also inserted into the
  # returned expression matrix.
  # ------------------------------------------------------------

  if (transform) {

    x_cluster <- asinh(x_cluster / cofactor)
    # also modify a whole matrix
    raw_ff[, channels] <- x_cluster
  }

  
  # ------------------------------------------------------------
  # 5. PCA
  #
  # PCA is performed on all selected fluorescence channels.
  # Centering and scaling give each channel equal variance before
  # dimensional reduction.
  # ------------------------------------------------------------

  pca_fit <- prcomp(
    x_cluster,
    center = TRUE,
    scale. = TRUE
  )
  
  xpca <- pca_fit$x
  
  
  # ------------------------------------------------------------
  # 6. Select number of PCs
  #
  # npcs = "all" retains all principal components.
  # A numeric value retains the first N PCs.
  # ------------------------------------------------------------

  if (!identical(npcs, "all")) {

    if (!is.numeric(npcs) || length(npcs) != 1 || npcs < 1) {
      stop(
        "`npcs` must be 'all' or a positive numeric value."
      )
    }
    
    npcs <- min(as.integer(npcs), ncol(xpca))
    xpca <- xpca[, seq_len(npcs), drop = FALSE]

  }
  

  # ------------------------------------------------------------
  # 7. UMAP
  #
  # UMAP is performed using the selected PCA coordinates rather
  # than the original fluorescence matrix.
  #
  # This makes the UMAP representation directly comparable to
  # downstream clustering performed in PCA space.
  # ------------------------------------------------------------
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  cat('X',nrow(x_cluster),"x",ncol(x_cluster))
  
  u_map <- uwot::umap(
    x_cluster,
    n_neighbors = n_neighbors,
    scale = umap_scale
  )
  
  colnames(u_map) <- c("UMAP1", "UMAP2")
  
  
  # ------------------------------------------------------------
  # 8. Return results
  # ------------------------------------------------------------
  cat(".\n")
  list(
    exprs = raw_ff,
    xpca = xpca,
    umap = u_map
  )
}


