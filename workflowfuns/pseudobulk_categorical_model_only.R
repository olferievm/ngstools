#' Test association between any categorical variable and pseudobulk gene expression
#'
#' Performs the analysis separately for each cell type using pseudobulk RNA-seq counts. Counts are filtered,
#' normalized with TMM, transformed using \code{voom}, and analyzed with
#' limma linear models.
#'
#' The categorical is treated as not ordered factorial variable.
#'
#' @param counts A matrix or data.frame of raw pseudobulk counts with
#'   samples in rows and genes in columns.
#' @param meta A data.frame containing one row per pseudobulk sample.
#'   Row names must match \code{counts}.
#' @param genes An optional data.frame with gene information and where all rownames equals colnames counts
#' @param var Character scalar giving the name of the studied categorical column
#'  values should factors.
#' @param contrasts Contrasts can be NULL or e.g. 'SLEvsHD'='DiagnosisSLE-DiagnosisHC'
#' @param cat_covariates Character vector of categorical covariates to
#'   include in the linear model if they contain more than one level.
#' @param num_covariates Character vector of categorical covariates to
#'   include in the linear model if they contain more than one level.
#' @param cell_types Optional character vector of cell types to analyze.
#'   Defaults to all cell types.
#' @param min.samples Minimum number of pseudobulk samples required for a
#'   cell type to be analyzed.
#'
#' @return
#' A data.frame containing differential expression statistics for all
#' analyzed cell types, including gene name, cell type, log fold change,
#' moderated t statistic, raw p-value, and FDR-adjusted p-value.
#'
#' @details
#' For each cell type the function:
#' \enumerate{
#'   \item Removes samples with missing genotype.
#'   \item Filters lowly expressed genes using edgeR::filterByExpr().
#'   \item Applies TMM normalization.
#'   \item Computes voom precision weights.
#'   \item Fits a limma linear model using an additive genotype effect.
#'   \item Returns differential expression statistics for the genotype coefficient.
#' }
#'
#' @importFrom edgeR DGEList filterByExpr calcNormFactors
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr bind_rows mutate select
#'
#' @export
#' 

pseudobulk_categorical_model_only <- function(counts, meta, genes=NULL, var, contrasts=NULL,
                                         cat_covariates = NULL, num_covariates = NULL,
                                         cell_column = "cell_type", cell_types = NULL,
                                         min_samples = 50, 
                                         p_value_cutoff =1, logFC_cutoff=0,
                                         min.count = 10, min.total.count = 100){
  
  stopifnot(is.matrix(counts) || is.data.frame(counts))
  stopifnot(is.data.frame(meta))
  stopifnot(var %in% colnames(meta))
  stopifnot(cell_column %in% colnames(meta))
  stopifnot(is.factor(meta[,var])|is.integer(meta[,var]))
  stopifnot(is.null(genes) || is.data.frame(genes))
  # Check cat_covariates are factor)
  if(!is.null(cat_covariates)){
    stopifnot(
      "All cat_covariates must be factors" = all(sapply(meta[, cat_covariates, drop = FALSE], is.factor))
    )
  }
  
  if(!is.null(num_covariates)){
    stopifnot(
      "All num_covariates must be numeric" = all(sapply(meta[, num_covariates, drop = FALSE], is.numeric))
    )
  }
  
  if (nrow(counts) != nrow(meta))
    stop("counts and meta contain different numbers of samples.")
  
  if (!identical(rownames(counts), rownames(meta)))
    stop("Row names of counts must match row names of meta.")
  
  if (!is.null(genes)){
    if(!identical(colnames(counts), rownames(genes)))
      stop("Col names of counts must match row names of genes.")
  }
  
  if(is.null(cell_types)){
    cell_types <- unique(meta$cell_type)
  }else{
    cell_types <- intersect(cell_types,unique(meta$cell_type))
  }
  
  screen_cell_types <- lapply(cell_types, function(cell_type){
    cat("Cells: ",cell_type,";", sep="")
    
    keep <- which(meta$cell_type == cell_type)
    cat(' samples:', length(keep))
    
    if(length(keep) < min_samples){cat("Not enough samples\n");return(NULL)}
    
    counts.ct <- counts[keep, ]
    meta.ct <- meta[keep, ]
    
    if(!is.null(genes)){
      dge <- DGEList(counts = t(counts.ct),
                     samples = meta.ct,
                     genes = genes)
    }else{
      dge <- DGEList(counts = t(counts.ct),
                     samples = meta.ct)
    }
    
    # remove missing values
    dge <- dge[,which(!is.na(dge$samples[,var]))]
    # remove empty rows
    dge <- dge[,which(dge$samples[,var] != "")]
    
    cat(' study var: ',var,', ', sep ="")
    
    dge$samples[,var] <- droplevels(factor(dge$samples[,var]))
    
    if (nlevels(dge$samples[,var]) < 2) {
      message(cell_type, ": only one level is present.\n")
      return(NULL)
    }
    
    gc <- table(dge$samples[,var])
    cat(paste(names(gc),gc, sep="="),sep=", ")
    
    f <- paste("~ 0 + ",var)
    
    # Add categorical covariate variables.
    if(!is.null(cat_covariates)){
      for(p in cat_covariates){
        f <- ifelse(levels(dge$samples[[p]]) > 1, paste(f, "+", p), f)
      }}
    
    # Add numerical covariates variables.
    if(!is.null(num_covariates)){
      for(p in num_covariates){
        f <- ifelse(length(unique(dge$samples[[p]])) > 1, paste(f, "+", p), f)
      }}
    
    cat(" formula:",f)
    
    keep.genes <- filterByExpr(dge,
                               design = model.matrix(as.formula(f),
                                                     dge$samples,
                                                     min.count = min.count,
                                                     min.total.count = min.total.count))
    
    cat(' genes: ', sum(keep.genes))
    
    dge <- dge[keep.genes,,keep.lib.sizes=FALSE]
    
    dge <- normLibSizes(dge, method="TMM")
    
    dge$samples <- droplevels(dge$samples)
    
    design <- model.matrix(as.formula(f), data= dge$samples)
    
    v <- voom(dge, design, plot=FALSE)
    
    fit <- lmFit(v, design)
      
    if (is.null(contrasts)) {
      
      coef.names <- colnames(design)
      
      coef.names <- coef.names[
        grepl(paste0("^", var), coef.names)
      ]
    
      fit <- eBayes(fit)
      
      return(list(voom = v, fit = fit))
      
    } else {
      
      contrast.matrix <- makeContrasts(
        contrasts = contrasts,
        levels = design
      )
      
      fit <- contrasts.fit(fit, contrast.matrix)
      fit <- eBayes(fit)
      return(list(voom = v, fit = fit))
      
    }
    
    cat("...Done\n")

    
  })
  
}
