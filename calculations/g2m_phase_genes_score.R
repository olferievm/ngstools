#' Calculate score for genes know to be regulated by cell cycle G2M phase
#'
#'
#' @param x A matrix of CPM values with gene names as rownames, or an edgeR::DGEList object.
#' @param controls A vector of column indices (or logicals) identifying control samples.
#'
#' @return A numeric vector of G2M phase gene response scores, one per sample.
#' @import edgeR
#' @export
#' 
#' @note 
#' 
#' @references

g2m_phase_genes_score <- function(x, controls = 1) {
  requireNamespace("edgeR", quietly = TRUE)
  
  # See above how the list was generated
  selgenes <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
                "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2",
                "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B",
                "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2",
                "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN",
                "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA")
  
  if (inherits(x, "DGEList")) {
    if (!"gene_name" %in% colnames(x$genes)) {
      stop("DGEList must have a 'gene_name' column in x$genes.")
    }
    
    matched_genes <- intersect(x$genes$gene_name, selgenes)
    
    if (length(matched_genes) == 0) stop("No response genes found. Abort.")
    if (length(matched_genes) < 10) stop("Only ", length(matched_genes), " genes found. Abort.")
    if (length(matched_genes) < 54) warning("NOTE: Found ", length(matched_genes), " genes.")
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
    if (length(matched_genes) < 54) warning("NOTE: Found ", length(matched_genes), " genes.")
    
    cpm_mat <- x[matched_genes, , drop = FALSE]
  } else {
    stop("Unsupported input type. Provide either a CPM matrix or a DGEList object.")
  }
  
  # Calculate control statistics
  # Need 2 genes
  if(nrow(cpm_mat) > 1){
    m_scaled <- scale(t(cpm_mat))
    m_score <- rowMeans(m_scaled)
    #hist(m_score, breaks = 100)
  }
  
  return(as.numeric(m_score))
}
