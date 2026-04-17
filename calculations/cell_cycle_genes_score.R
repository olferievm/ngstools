#' Calculate score for genes know to be regulated by cell cycle S phase
#'
#'
#' @param x A matrix of CPM values with gene names as rownames, or an edgeR::DGEList object.
#' @param ref.genes Select either "PBMC" manually selected based on PBMC data or "Seurat" provided for scRNA-seq
#' @return A numeric vector of S phase gene response scores, one per sample.
#' @import edgeR
#' @export
#' 
#' @note 
#' 
#' @references

cell_cycle_genes_score <- function(x, ref.genes = c("PBMC","Suerat")) {
  requireNamespace("edgeR", quietly = TRUE)

  if(ref.genes[1] == "Seurat"){
  # See above how the list was generated. 
  s_genes <- c(
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
    "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS", "RFC2", "RPA2", "NASP",
    "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2",
    "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1",
    "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"
  )
  g2m_genes <- c(
  "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
  "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2",
  "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B",
  "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2",
  "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN",
  "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA")
  }else{
  # Selected on PCA analysis in PBMC
  s_genes <- c("HELLS", "MSH2", "PCNA", "FEN1", "WDR76", "PRIM1",
               "MCM4", "CHAF1B", "RRM1", "MCM6", "POLA1", "MCM2", "UNG", "MCM5")
  g2m_genes <- c("MKI67", "TOP2A", "CDK1", "CCNB2", "CDC20", "CDC25C",
                   "BUB1", "BUB1B", "AURKA", "AURKB", "KIF11", "KIF20B",
                   "KIF23", "TPX2", "NUSAP1", "CENPA", "CENPE",
                   "HJURP", "DLGAP5", "ANLN")
  }
  
  
  selgenes <- c(s_genes, g2m_genes)
  
  if (inherits(x, "DGEList")) {
    if (!"gene_name" %in% colnames(x$genes)) {
      stop("DGEList must have a 'gene_name' column in x$genes.")
    }
    
    matched_genes <- intersect(x$genes$gene_name, selgenes)
 
    if (length(matched_genes) == 0) stop("No response genes found. Abort.")
    if (length(matched_genes) < 10) stop("Only ", length(matched_genes), " genes found. Abort.")
    if (length(matched_genes) < 45) warning("NOTE: Found ", length(matched_genes), " genes.")
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
    if (length(matched_genes) < 45) warning("NOTE: Found ", length(matched_genes), " genes.")
    
    cpm_mat <- x[matched_genes, , drop = FALSE]
  } else {
    stop("Unsupported input type. Provide either a CPM matrix or a DGEList object.")
  }
  
  
  s_genes <- intersect(rownames(cpm_mat), s_genes)
  g2m_genes <- intersect(rownames(cpm_mat), g2m_genes)
  
  # Need at least 5 genes
  if(length(s_genes) < 5) stop("Few genes, impossible to estimate. Abort.")
  if(length(g2m_genes) < 5) stop("Few genes, impossible to estimate. Abort.")
  
    pca <- FactoMineR::PCA(t(cpm_mat), graph = FALSE)
    pc_loads <- pca$var$coord[,1:2]
    
    score_set <- function(gene_set, comp) {
      mean(abs(pc_loads[intersect(rownames(pc_loads), gene_set), comp]))
    }
    
    assoc <- data.frame(
      Component = c("Dim.1","Dim.2"),
      S = c(score_set(s_genes2,"Dim.1"), score_set(s_genes2,"Dim.2")),
      G2M = c(score_set(g2m_genes2,"Dim.1"), score_set(g2m_genes2,"Dim.2"))
    )
    
    S_comp  <- assoc$Component[which.max(assoc$S)]
    G2M_comp <- assoc$Component[which.min(assoc$S)]
    
    sample_scores <- pca$ind$coord
    
    S_score  <- sample_scores[, S_comp]
    G2M_score <- sample_scores[, G2M_comp]
    
    # tmp <- data.frame(S_score = S_score, G2M_score = G2M_score)
    # g <- ggplot2::ggplot(tmp, aes(x = S_score, y = G2M_score)) + geom_point()
    # g
    #hist(G2M_score, breaks = 100)
    return(data.frame(s_phase_score = S_score, g2m_phase_score = G2M_score))

}