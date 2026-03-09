#' Calculate score for genes know to be regulated by TNF
#'
#'
#' @param x A matrix of CPM values with gene names as rownames, or an edgeR::DGEList object.
#' @param controls A vector of column indices (or logicals) identifying control samples.
#'
#' @return A numeric vector of TNF gene response scores, one per sample.
#' @import edgeR
#' @export
#' 
#' @note The list of TNF induced genes was generated as follow: 1. STRING interaction network of TNF co-expressed genes, 2. GENEMANIA interaction network of coexpressed genes. 3. GSE89408 synovial biopsy dataset of OA and RA patients known to have genes induced by chronic TNF signaling. 4.ChatGPT request on top correlating genes. Strong / direct TNF-responsive or TNF-pathway genes
#' @note NFKBIE (IκBε) Classic NF-κB feedback inhibitor. Rapidly induced by TNF to dampen NF-κB signaling.
#' @note TNFAIP3 (A20) One of the strongest TNF response genes. Induced within minutes; terminates TNFR signaling via ubiquitin editing. Central in RA/OA inflammation.
#' @note BIRC3 (cIAP2) Direct TNF-α target. Anti-apoptotic, TNFR1-proximal signaling component; frequently upregulated in TNF-stimulated synovial fibroblasts.
#' @note DUSP2 Induced by TNF; negative regulator of MAPK signaling downstream of TNFR. Highly expressed in activated T cells and inflamed synovium.
#' @note IER5 Immediate early gene induced by TNF/NF-κB and stress signaling; less specific but robust.
#' @note SLC2A6 (GLUT6) TNF-induced metabolic reprogramming gene in macrophages and fibroblasts; commonly elevated in inflammatory synovium.
#' @note LTB (Lymphotoxin-β) Member of the TNF superfamily. Co-regulated with TNF signaling; reflects sustained inflammatory circuits.
#' @note TAGAP Upregulated in TNF-driven T-cell activation; strong autoimmune disease association (RA, SLE).
#' @note These are the most defensible if you need mechanistic credibility:
#' @note TNFAIP8L2 (TIPE2) Direct TNF-α–induced gene; negative regulator of TNF/NF-κB signaling, macrophage activation, inflammation.
#' @note IRF8 TNF activates IRF8 in myeloid cells; critical for TNF-driven inflammatory and IFN-associated transcriptional programs.
#' @note SP140 Induced downstream of TNF/NF-κB in immune cells; enriched in inflammatory and autoimmune contexts (including RA-like biology).
#' @note PSME2 (PA28β) - TNF-α induces immunoproteasome components; part of TNF/IFN-driven antigen-processing response.
#' @note KMO - Upregulated by TNF in macrophages and synovial fibroblasts; links TNF signaling to kynurenine metabolism and inflammation.
#' @note COTL1 TNF-α induces actin remodeling genes; COTL1 increases in TNF-stimulated macrophages and synovial cells.
#' @note LPXN TNF-responsive focal adhesion protein; linked to TNF-driven migration and adhesion in inflammatory cells.
#' @note PFN1 Cytoskeletal remodeling downstream of TNF/NF-κB signaling.
#' @note RAB8A TNF affects vesicular trafficking; RAB8A participates in TNF-regulated membrane dynamics.
#' 
#' @references

TNF_genes_score <- function(x, controls = 1) {
  requireNamespace("edgeR", quietly = TRUE)
  
  # See above how the list was generated
  tnf_genes <- c(
    "NFKBIE", "TNFAIP3", "BIRC3", "SLC2A6", "DUSP2", "IER5", "TAGAP", "LTB",
    "TNFAIP8L2", "IRF8", "SP140", "PSME2", "KMO", "COTL1", "LPXN", "PTPN6"
  )
  
  if (inherits(x, "DGEList")) {
    if (!"gene_name" %in% colnames(x$genes)) {
      stop("DGEList must have a 'gene_name' column in x$genes.")
    }
    
    matched_genes <- intersect(x$genes$gene_name, tnf_genes)
    
    if (length(matched_genes) == 0) stop("No TNF response genes found. Abort.")
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
    matched_genes <- intersect(rownames(x), tnf_genes)
    
    if (length(matched_genes) == 0) stop("No TNF response genes found. Abort.")
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
