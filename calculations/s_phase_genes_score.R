#' Calculate score for genes know to be regulated by cell cycle S phase
#'
#'
#' @param x A matrix of CPM values with gene names as rownames, or an edgeR::DGEList object.
#' @param controls A vector of column indices (or logicals) identifying control samples.
#'
#' @return A numeric vector of S phase gene response scores, one per sample.
#' @import edgeR
#' @export
#' 
#' @note 
#' 
#' @references

s_phase_genes_score <- function(x, controls = 1) {
  requireNamespace("edgeR", quietly = TRUE)
  
  # See above how the list was generated
  selgenes <- c(
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"
  )
  
  if (inherits(x, "DGEList")) {
    if (!"gene_name" %in% colnames(x$genes)) {
      stop("DGEList must have a 'gene_name' column in x$genes.")
    }
    
    matched_genes <- intersect(x$genes$gene_name, selgenes)
    
    if (length(matched_genes) == 0) stop("No response genes found. Abort.")
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
    matched_genes <- intersect(rownames(x), selgenes)
    
    if (length(matched_genes) == 0) stop("No genes found. Abort.")
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
