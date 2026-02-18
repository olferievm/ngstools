#' Calculate sex specific gene score from expression matrix or edgeR object
#'
#'
#' @param x A matrix of CPM values with gene names as rownames, or an edgeR::DGEList object.
#'
#' @return A numeric vector of interferon scores, one per sample.
#' @import edgeR
#' @export
#' 
#' 
#' @references
#'

sexspecificgenes <- function(x) {
  requireNamespace("edgeR", quietly = TRUE)

  m_score <- NULL
  f_score <- NULL
  sex_score <- NULL
  
  # Define the 21 interferon-stimulated genes
  male_genes <- c(
       "TMSB4Y", "EIF1AY", "UTY", "NLGN4Y", "KDM5D", "DDX3Y", "TXLNGY",
       "RPS4Y1", "RPS4Y2", "USP9Y", "PRKY", "TTTY15", "TTTY14", "TTTY10",
       "LINC00278", "ZFY", "LOC107987338", "LOC102724150", "LOC105377223",
       "SRY", "LOC107987348", "LOC107987350", "BCORP1", "HSFY1", "HSFY2",
       "FAM224A", "FAM224B", "LOC105377225", "TBL1Y")
  
  female_genes <- c("XIST", "TSIX", "KDM6A", "ZFX", "PUDP", "PRKX",
                    "DDX3X", "JPX", "TXLNG", "SMC1A", "EIF1AX",
                    "RPS4X", "ZRSR2", "EIF2S3", "ARSD",
                    "ALG13", "TRAPPC2", "SYAP1", "KDM5C",
                    "CXorf38", "USP9X")

  sex_genes <- c(female_genes, male_genes)
  
  if (inherits(x, "DGEList")) {
    if (!"gene_name" %in% colnames(x$genes)) {
      stop("DGEList must have a 'gene_name' column in x$genes.")
    }
    
    matched_genes <- intersect(x$genes$gene_name, sex_genes)
    
    if (length(matched_genes) == 0) stop("No sex specific genes found. Abort.")
    if (length(matched_genes) < 10) stop("Only ", length(matched_genes), " genes found. Abort.")
    if (length(matched_genes) < 21) warning("NOTE: Found ", length(matched_genes), " genes.")
    if (anyDuplicated(matched_genes)) stop("Duplicated gene names found. Abort.")
    
    gene_indices <- which(x$genes$gene_name %in% matched_genes)
    gene_labels <- x$genes$gene_name[gene_indices]
    
    dup_names <- duplicated(gene_labels)
    
    gene_indices <- gene_indices[!dup_names]
    gene_labels <- gene_labels[!dup_names] 
    
    
    cpm_mat <- edgeR::cpm(x, normalized.lib.sizes = TRUE,
                          log = TRUE,
                          prior.count = 2)
    
    cpm_mat <- cpm_mat[gene_indices, , drop = FALSE]
    rownames(cpm_mat) <- gene_labels
    
  } else if (is.matrix(x)) {
    matched_genes <- intersect(rownames(x), interferon_genes)
    
    if (length(matched_genes) == 0) stop("No sex specific genes found. Abort.")
    if (length(matched_genes) < 10) stop("Only ", length(matched_genes), " genes found. Abort.")
    if (length(matched_genes) < 21) warning("NOTE: Found ", length(matched_genes), " genes.")
    
    cpm_mat <- x[matched_genes, , drop = FALSE]
    
  } else {
    stop("Unsupported input type. Provide either a CPM matrix or a DGEList object.")
  }
  
  
  male_genes <- intersect(rownames(cpm_mat), male_genes)
  female_genes <- intersect(rownames(cpm_mat), female_genes)
  

  
  # Need 2 genes
  if(length(male_genes) > 1){
    m_cpm <- cpm_mat[male_genes, ,drop = FALSE]
    m_scaled <- scale(t(m_cpm))
    m_score <- rowMeans(m_scaled)
    #hist(m_score, breaks = 100)
  }
 # Need 2 genes
  if(length(female_genes) > 1){
    f_cpm <- cpm_mat[female_genes, ,drop = FALSE]
    f_scaled <- scale(t(f_cpm))
    f_score <- rowMeans(f_scaled)
    #hist(f_score, breaks = 100)
  }
  
  
  if(length(male_genes) > 1 & length(female_genes) > 1){
    f_scaled <- f_scaled * (-1)
    a_scaled <- cbind(m_scaled,f_scaled)
    sex_score <- rowMeans(a_scaled)
    #hist(sex_score, breaks = 100)
  }
  # Calculate final score statistics
   
   return(
     list(
       male_gene_score = m_score,
       female_gene_score = f_score,
       sex_gene_score = sex_score
     )
   )
}
