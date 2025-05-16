#' Calculate interferon score from expression matrix or edgeR object
#'
#'
#'
#'
#'
#' @param x A matrix of CPM values with gene names as rownames, or an edgeR::DGEList object.
#' @param controls A vector of column indices (or logicals) identifying control samples.
#'
#' @return A numeric vector of interferon scores, one per sample.
#' @import edgeR
#' @export


interferon_21genes <- function(x, controls = 1) {
  requireNamespace("edgeR", quietly = TRUE)
  
  # Define the 21 interferon-stimulated genes
  interferon_genes <- c(
    "SPATS2L", "ISG15", "IFIT1", "IFIT2", "IFIT3",
    "RSAD2", "HERC5", "IFI44", "CMPK2", "IFI6", 
    "OAS1", "OAS2", "OASL", "OAS3", "MX1",
    "EPSTI1", "LY6E", "HERC6", "LAMP3", "IFI27", "OTOF"
  )
  
  if (inherits(x, "DGEList")) {
    if (!"gene_name" %in% colnames(x$genes)) {
      stop("DGEList must have a 'gene_name' column in x$genes.")
    }
    
    matched_genes <- intersect(x$genes$gene_name, interferon_genes)
    
    if (length(matched_genes) == 0) stop("No interferon genes found. Abort.")
    if (length(matched_genes) < 10) stop("Only ", length(matched_genes), " genes found. Abort.")
    if (length(matched_genes) < 21) warning("NOTE: Found ", length(matched_genes), " genes.")
    if (anyDuplicated(matched_genes)) stop("Duplicated gene names found. Abort.")
    
    gene_indices <- which(x$genes$gene_name %in% matched_genes)
    gene_labels <- x$genes$gene_name[gene_indices]
    
    cpm_mat <- edgeR::cpm(x, normalized.lib.sizes = TRUE,
                          log = TRUE,
                          prior.count = 2)
    
    cpm_mat <- cpm_mat[gene_indices, , drop = FALSE]
    rownames(cpm_mat) <- gene_labels
    
  } else if (is.matrix(x)) {
    matched_genes <- intersect(rownames(x), interferon_genes)
    
    if (length(matched_genes) == 0) stop("No interferon genes found. Abort.")
    if (length(matched_genes) < 10) stop("Only ", length(matched_genes), " genes found. Abort.")
    if (length(matched_genes) < 21) warning("NOTE: Found ", length(matched_genes), " genes.")
    
    cpm_mat <- x[matched_genes, , drop = FALSE]
  } else {
    stop("Unsupported input type. Provide either a CPM matrix or a DGEList object.")
  }
  
  # Calculate control statistics
  control_means <- rowMeans(cpm_mat[, controls, drop = FALSE])
  control_sds   <- apply(cpm_mat[, controls, drop = FALSE], 1, sd)
  
  # Standardize expression: (x - mean) / sd
  z_scores <- sweep(cpm_mat, 1, control_means, "-")
  z_scores <- sweep(z_scores, 1, control_sds, "/")
  
  # Interferon score is the average Z-score per sample
  colMeans(z_scores, na.rm = TRUE)
}
