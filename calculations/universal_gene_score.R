#' Calculate gene score from expression matrix or edgeR object.
#'
#'
#' @param x A matrix of CPM values with gene names as rownames, or an edgeR::DGEList object.
#' @param genes A vector of genes used for generating a score, e.g. genes of interferon response.
#' @param controls A vector of column indices (or logicals) identifying control samples.
#'
#' @return A numeric vector of interferon scores, one per sample.
#' @import edgeR
#' @export
#' 
#' 
#' @references Yao, Y., Higgs, B. W., Morehouse, C., de los Reyes, M., Trigona, W., Brohawn, P., White, W., Zhang, J., White, B., Coyle, A. J., Kiener, P. A., & Jallal, B. (2009). Development of Potential Pharmacodynamic and Diagnostic Markers for Anti-IFN-α Monoclonal Antibody Trials in Systemic Lupus Erythematosus. Human Genomics and Proteomics, 1(1). https://doi.org/10.4061/2009/374312
#'

universal_gene_score <- function(x, genes, method = c("average", "pca","z-score"), controls = 1) {
  requireNamespace("edgeR", quietly = TRUE)
  
  if (inherits(x, "DGEList")) {
    if (!"gene_name" %in% colnames(x$genes)) {
      stop("DGEList must have a 'gene_name' column in x$genes.")
    }
    
    matched_genes <- intersect(x$genes$gene_name, genes)
    
    if (length(matched_genes) == 0) stop("No selected genes found. Abort.")
    if (anyDuplicated(matched_genes)) stop("Duplicated gene names found. Abort.")
    if (length(matched_genes) > 0) cat("NOTE: Found ", length(matched_genes), " genes. ")
    
    gene_indices <- which(x$genes$gene_name %in% matched_genes)
    gene_labels <- x$genes$gene_name[gene_indices]
    
    cpm_mat <- edgeR::cpm(x, normalized.lib.sizes = TRUE,
                          log = TRUE,
                          prior.count = 2)
    
    cpm_mat <- cpm_mat[gene_indices, , drop = FALSE]
    rownames(cpm_mat) <- gene_labels
    
  } else if (is.matrix(x)) {
    
    matched_genes <- intersect(rownames(x), interferon_genes)
    
    if (length(matched_genes) == 0) stop("No selected genes found. Abort.")
    if (anyDuplicated(matched_genes)) stop("Duplicated gene names found. Abort.")
    if (length(matched_genes) > 0) cat("NOTE: Found ", length(matched_genes), " genes. ")
    
    cpm_mat <- x[matched_genes, , drop = FALSE]
  } else {
    stop("Unsupported input type. Provide either a CPM matrix or a DGEList object.")
  }
  
  if(method == 'average'){
    z_scores <- scale(t(cpm_mat), center = TRUE, scale = TRUE)
    return(rowMeans(z_scores))
  }
  
  if(method == 'z-score'){
  # Calculate control statistics
     control_means <- rowMeans(cpm_mat[, controls, drop = FALSE])
     control_sds   <- apply(cpm_mat[, controls, drop = FALSE], 1, sd)
  
  # Standardize expression: (x - mean) / sd
     z_scores <- sweep(cpm_mat, 1, control_means, "-")
     z_scores <- sweep(z_scores, 1, control_sds, "/")
  
  # Interferon score is the average Z-score per sample
     return(colMeans(z_scores, na.rm = TRUE))
  }
  
  if(method == 'pca'){
    pca <- prcomp(t(cpm_mat), center = TRUE, scale. = TRUE)
    eigengene <- pca$x[, 1]
#   The sign of PC1 is arbitrary.
#   To make higher eigengene values correspond to higher average expression.
    
    if(cor(
      eigengene,
      colMeans(scale(t(sub)))
    ) < 0)
    {
      eigengene <- -eigengene
    }
    return(eigengene)
  }
  
}
