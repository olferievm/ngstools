tpm <- function(y, log = TRUE, prior = 2) {
  # y: a DGEList object from edgeR
  # log: whether to log2-transform the TPM output
  # prior: pseudocount added before log transformation (to avoid log(0))
  
  # Validate input is a DGEList
  if (!inherits(y, "DGEList")) {
    stop("Provide a DGEList object.")
  }
  
  # Check gene lengths exist and match number of rows
  if (is.null(y$genes$length) || length(y$genes$length) != nrow(y)) {
    stop("DGEList object must include a 'genes$length' column matching the number of genes.")
  }
  
  # Compute Reads Per Kilobase (RPK)
  rpk <- y$counts / y$genes$length
  
  # Compute per-sample scaling factors (sum of RPKs)
  sf <- colSums(rpk)
  
  # Calculate TPM
  tpm <- t(t(rpk) / sf) * 1e6
  
  # Return either raw TPM or log2-transformed TPM
  if (!log) {
    return(tpm)
  } else {
    return(log2(tpm + prior))
  }
}
