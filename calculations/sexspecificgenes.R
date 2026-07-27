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
  
  # Define the 21 interferon-stimulated genes
  male_genes <- c(
       "TMSB4Y", "EIF1AY", "UTY", "NLGN4Y", "KDM5D", "DDX3Y", "TXLNGY",
       "RPS4Y1", "RPS4Y2", "USP9Y", "PRKY", "TTTY15", "TTTY14", "TTTY10",
       "LINC00278", "ZFY", "LOC107987338", "LOC102724150", "LOC105377223",
       "SRY", "LOC107987348", "LOC107987350", "BCORP1", "HSFY1", "HSFY2",
       "FAM224A", "FAM224B", "LOC105377225", "TBL1Y") %>% unique()
  
  female_genes <- c("XIST", "TSIX", "KDM6A", "ZFX", "PRKX",
                    "DDX3X", "JPX", "TXLNG", "EIF1AX",
                    "RPS4X", "USP9X") %>% unique()

  sex_genes <- c(female_genes, male_genes)
  
  if (inherits(x, "DGEList")) {
    if (!"gene_name" %in% colnames(x$genes)) {
      stop("DGEList must have a 'gene_name' column in x$genes.")
    }
    
        if(anyDuplicated(x$genes$gene_name)){
            cat("Duplicated gene names found")
            x <- x[!duplicated(x$genes$gene_name),]
         }
    
       matched_male_genes <- intersect(x$genes$gene_name, male_genes)
       matched_female_genes <- intersect(x$genes$gene_name, female_genes)
       matched_genes <- intersect(x$genes$gene_name, sex_genes)
    

    if(length(matched_male_genes) < 2 & length(matched_female_genes < 2)){cat("No male of female specific genes found\n"); return(NULL)}
    if(length(matched_male_genes) < 2){cat("No male specific genes found\n")}
    if(length(matched_female_genes) < 2){cat("No female specific genes found\n")}

    cpm_mat <- edgeR::cpm(x, normalized.lib.sizes = TRUE,
                          log = TRUE,
                          prior.count = 2)
    
    gene_indices <- which(x$genes$gene_name %in% matched_genes)
    gene_labels <- x$genes$gene_name[gene_indices]
    cpm_mat <- cpm_mat[gene_indices, , drop = FALSE]
    rownames(cpm_mat) <- gene_labels
    
  } else if (is.matrix(x)) {
    
    if(anyDuplicated(rownames(x))){
      cat("Duplicated gene names found")
      x <- x[,!duplicated(rownames(x))]
    }
    
    genes <- rownames(x)
    matched_male_genes <- intersect(genes, male_genes)
    matched_female_genes <- intersect(genes, female_genes)
    matched_genes <- intersect(genes, sex_genes)
    
    if(length(matched_male_genes) < 2 & length(matched_female_genes < 2)){cat("No male of female specific genes found\n"); return(NULL)}
    if(length(matched_male_genes) < 2){cat("No male specific genes found\n")}
    if(length(matched_female_genes) < 2){cat("No female specific genes found\n")}
    
    gene_indices <- which(x$genes$gene_name %in% matched_genes)
    gene_labels <- x$genes$gene_name[gene_indices]
    cpm_mat <- x[matched_genes, , drop = FALSE]
    
  } else {
    stop("Unsupported input type. Provide either a CPM matrix or a DGEList object.")
  }
  
  # Need 2 genes
  if(length(matched_male_genes) > 1){
    m_cpm <- cpm_mat[matched_male_genes, ,drop = FALSE]
    m_scaled <- scale(t(m_cpm))
    m_score <- rowMeans(m_scaled)
    #hist(m_score, breaks = 100)
  }
 # Need 2 genes
  if(length(matched_female_genes) > 1){
    f_cpm <- cpm_mat[matched_female_genes, ,drop = FALSE]
    f_scaled <- scale(t(f_cpm))
    f_score <- rowMeans(f_scaled)
    #hist(f_score, breaks = 100)
  }
  
  
  if(length(matched_male_genes) > 1 & length(matched_female_genes) > 1){
    f_scaled <- f_scaled * (-1)
    a_scaled <- cbind(m_scaled, f_scaled)
    sex_score <- rowMeans(a_scaled)
    #hist(sex_score, breaks = 100)
  }
  # Calculate final score statistics
   
   return(
     list(
       male_genes_score = m_score,
       female_genes_score = f_score,
       sex_genes_score = sex_score
     )
   )
}
