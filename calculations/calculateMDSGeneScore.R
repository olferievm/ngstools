#' Calculate Gene Score Matrix (e.g. correlation with MDS axes or SVs)
#'
#' @param y edgeR::DGEList object. Must include `genes` slot with gene annotations.
#' @param ma Matrix (samples x components) to correlate with (e.g., mds$eigenvectors or sva$sv)
#' @param top Integer. Number of top genes to select (used only for 'common' or 'pairwise' gene selection).
#' @param gene_selection One of 'none', 'common', 'pairwise'. Method for selecting top variable genes.
#' @param gene_subset Optional numeric or logical vector to subset rows of y before calculating scores.
#'
#' @return A tibble: gene annotation (from y$genes) joined with gene score matrix.
#' @example mds <- plotMDS(dge, returnData = TRUE)
#' @example gene_score_df <- calculateMDSGeneScore(dge, ma = mds$eigenvectors, top = 500, gene_selection = 'common')
#'
#'


calculateMDSGeneScore <- function(y, ma, top = 500, 
                               gene_selection = c('none', 'common', 'pairwise'), 
                               gene_subset = 1:nrow(y)) {
  # Check and set gene selection
  gene_selection <- match.arg(gene_selection)
  
  # Subset y if needed
  y <- y[gene_subset, , keep.lib.sizes = TRUE]
  
  # Correct rownames in case are empy (avoid duplicated names)
  if(any(duplicated(y$genes$gene))){stop("Duplicated y$genes$gene - are not allowed")}
  rownames(y$genes) <- y$genes$gene
  rownames(y$counts) <- y$genes$gene
  
  # Compute logCPM (genes x samples)
  logCPM <- edgeR::cpm(y, normalized.lib.sizes = TRUE, log = TRUE)
  
  # Initialize top_genes
  top_genes <- NULL
  
  # Select top variable genes (common)
  if (gene_selection == 'common') {
    gene_vars <- apply(logCPM, 1, var)
    top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:top]
  }
  
  # Select top genes by pairwise logFC variability
  if (gene_selection == 'pairwise') {
    lfc_sums <- numeric(nrow(logCPM))
    names(lfc_sums) <- rownames(logCPM)
    n_samples <- ncol(logCPM)
    
    for (i in 1:(n_samples - 1)) {
      for (j in (i + 1):n_samples) {
        lfc <- logCPM[, i] - logCPM[, j]
        lfc_sums <- lfc_sums + abs(lfc)
      }
    }
    
    top_genes <- names(sort(lfc_sums, decreasing = TRUE))[1:top]
  }
  
  # If no selection, use all genes
  if (gene_selection == 'none') {
    top_genes <- rownames(logCPM)
  }
  
  if (is.null(top_genes) || length(top_genes) < 3) {
    stop("No sufficient genes selected to calculate gene score.")
  }
  
  # Subset and transpose logCPM (samples x selected genes)
  logCPM_tr <- t(logCPM[top_genes, ])
  
  # Check matrix dimensions
  if (nrow(logCPM_tr) != nrow(ma)) {
    stop("Number of samples in CPM and score matrix (ma) do not match.")
  }
  
  # Compute correlations: genes (rows) × dimensions (columns)
  gene_scores <- cor(logCPM_tr, ma)
  
  # Format gene names
  gene_scores <- as.data.frame(gene_scores)
  gene_scores <- tibble::rownames_to_column(gene_scores, var = "gene")
  
  # Extract gene annotations from DGEList
  if (is.null(y$genes)) {
    stop("DGEList object must contain a `genes` annotation slot.")
  }
  
  gene_annots <- y$genes[top_genes, , drop = FALSE]
  
  # Join and return
  dplyr::left_join(gene_annots, gene_scores, by = "gene")
}





