#' Test gene-expression associations in pseudobulk RNA-seq
#'
#' Performs a genome-wide gene-expression association analysis separately
#' for each cell type using pseudobulk RNA-seq counts. Selected gene
#' expression values already present in `meta` are used as continuous
#' predictors. Counts are filtered, TMM normalized, transformed using
#' voom, and analyzed with limma linear models.
#'
#' For each predictor gene, the model is:
#'
#'   response_gene ~ predictor_gene + covariates
#'
#' where response_gene is the voom logCPM expression of every retained
#' gene in the pseudobulk count matrix.
#'
#' @param counts A matrix or data.frame of raw pseudobulk counts with
#'   samples in rows and genes in columns.
#' @param meta A data.frame containing one row per pseudobulk sample.
#'   Row names must match `counts`.
#' @param genes An optional data.frame with gene information. Row names
#'   must match colnames(counts).
#' @param predictor_genes Character vector giving names of columns in
#'   `meta` containing log-transformed expression values of genes to test.
#' @param cat_covariates Character vector of categorical covariates.
#' @param num_covariates Character vector of numerical covariates.
#' @param cell_column Name of the cell-type column in `meta`.
#' @param cell_types Optional character vector of cell types to analyze.
#' @param min_samples Minimum number of pseudobulk samples required.
#' @param p_value_cutoff P-value cutoff passed to topTable.
#' @param logFC_cutoff Minimum absolute association coefficient passed
#'   to topTable.
#' @param scale_predictors Logical. If TRUE, predictor expression is
#'   standardized within each cell type.
#'
#' @return
#' A data.frame containing association statistics for all predictor genes
#' and cell types.
#'
#' @export
#'

pseudobulk_gene_association <- function(
    counts,
    meta,
    genes = NULL,
    predictor_genes,
    cat_covariates = NULL,
    num_covariates = NULL,
    cell_column = "cell_type",
    cell_types = NULL,
    min_samples = 50,
    p_value_cutoff = 1,
    logFC_cutoff = 0,
    scale_predictors = TRUE,
    min.count = 10,
    min.total.count = 100
) {
  
  stopifnot(is.matrix(counts) || is.data.frame(counts))
  stopifnot(is.data.frame(meta))
  stopifnot(cell_column %in% colnames(meta))
  stopifnot(all(predictor_genes %in% colnames(meta)))
  
  if (!is.null(genes)) {
    stopifnot(is.data.frame(genes))
    
    if (!identical(colnames(counts), rownames(genes))) {
      stop("Col names of counts must match row names of genes.")
    }
  }
  
  if (!is.null(cat_covariates)) {
    stopifnot(
      "All cat_covariates must be factors" =
        all(sapply(
          meta[, cat_covariates, drop = FALSE],
          is.factor
        ))
    )
  }
  
  if (!is.null(num_covariates)) {
    stopifnot(
      "All num_covariates must be numeric" =
        all(sapply(
          meta[, num_covariates, drop = FALSE],
          is.numeric
        ))
    )
  }
  
  if (nrow(counts) != nrow(meta)) {
    stop("counts and meta contain different numbers of samples.")
  }
  
  if (!identical(rownames(counts), rownames(meta))) {
    stop("Row names of counts must match row names of meta.")
  }
  
  if (is.null(cell_types)) {
    cell_types <- unique(meta[[cell_column]])
  } else {
    cell_types <- intersect(
      cell_types,
      unique(meta[[cell_column]])
    )
  }
  
  results <- list()
  result_index <- 1
  
  for (cell_type in cell_types) {
    
    cat("\n", cell_type, ": ", sep = "")
    
    keep <- which(meta[[cell_column]] == cell_type)
    
    cat("samples:", length(keep))
    
    if (length(keep) < min_samples) {
      cat(" -- skipped\n")
      next
    }
    
    counts.ct <- counts[keep, , drop = FALSE]
    meta.ct <- meta[keep, , drop = FALSE]
    
    #
    # Process each predictor gene separately.
    #
    
    for (predictor in predictor_genes) {
      
      cat("\n  ", predictor, ": ", sep = "")
      
      #
      # Remove samples with missing predictor.
      #
      
      keep.predictor <- !is.na(meta.ct[[predictor]])
      
      counts.pred <- counts.ct[keep.predictor, , drop = FALSE]
      meta.pred <- meta.ct[keep.predictor, , drop = FALSE]
      
      if (nrow(meta.pred) < min_samples) {
        cat("too few samples")
        next
      }
      
      #
      # Check that predictor has variation.
      #
      
      if (length(unique(meta.pred[[predictor]])) < 2) {
        cat("no predictor variation")
        next
      }
      
      #
      # Standardize predictor if requested.
      #
      # This makes the coefficient represent the change in response
      # associated with a one-SD increase in predictor expression.
      #
      
      if (scale_predictors) {
        meta.pred$predictor_expression <-
          as.numeric(scale(meta.pred[[predictor]]))
      } else {
        meta.pred$predictor_expression <-
          meta.pred[[predictor]]
      }
      
      #
      # Build initial design.
      #
      
      formula_terms <- "predictor_expression"
      
      #
      # Add categorical covariates.
      #
      
      if (!is.null(cat_covariates)) {
        
        for (p in cat_covariates) {
          
          if (nlevels(droplevels(meta.pred[[p]])) > 1) {
            
            formula_terms <- paste(
              formula_terms,
              "+",
              p
            )
          }
        }
      }
      
      #
      # Add numerical covariates.
      #
      
      if (!is.null(num_covariates)) {
        
        for (p in num_covariates) {
          
          if (length(unique(
            meta.pred[[p]][!is.na(meta.pred[[p]])]
          )) > 1) {
            
            formula_terms <- paste(
              formula_terms,
              "+",
              p
            )
          }
        }
      }
      
      formula <- as.formula(
        paste("~", formula_terms)
      )
      
      design <- model.matrix(
        formula,
        data = meta.pred
      )
      
      #
      # Construct DGEList.
      #
      
      if (!is.null(genes)) {
        
        dge <- edgeR::DGEList(
          counts = t(counts.pred),
          samples = meta.pred,
          genes = genes
        )
        
      } else {
        
        dge <- edgeR::DGEList(
          counts = t(counts.pred),
          samples = meta.pred
        )
      }
      
      #
      # Filter lowly expressed genes.
      #
      
      keep.genes <- edgeR::filterByExpr(
        dge,
        design = design,
        min.count = min.count,
        min.total.count = min.total.count
      )
      
      cat(
        "genes:",
        sum(keep.genes),
        " "
      )
      
      dge <- dge[
        keep.genes,
        ,
        keep.lib.sizes = FALSE
      ]
      
      #
      # TMM normalization.
      #
      
      dge <- edgeR::normLibSizes(
        dge,
        method = "TMM"
      )
      
      #
      # Voom transformation.
      #
      
      v <- limma::voom(
        dge,
        design,
        plot = FALSE
      )
      
      #
      # Fit limma model.
      #
      
      fit <- limma::lmFit(
        v,
        design
      )
      
      fit <- limma::eBayes(
        fit
      )
      
      #
      # Extract association statistics.
      #
      
      res <- limma::topTable(
        fit,
        coef = "predictor_expression",
        number = Inf,
        p.value = p_value_cutoff,
        lfc = logFC_cutoff,
        sort.by = "P"
      ) |>
        as.data.frame()
      
      #
      # Add gene name if not already present.
      #
      
      if (!"gene_name" %in% colnames(res)) {
        
        res <- tibble::rownames_to_column(
          res,
          var = "gene_name"
        )
      }
      
      #
      # Add predictor and cell type.
      #
      
      res <- res |>
        dplyr::mutate(
          cell_type = cell_type,
          predictor = predictor
        ) |>
        dplyr::select(
          predictor,
          cell_type,
          gene_name,
          dplyr::everything()
        )
      
      results[[result_index]] <- res
      result_index <- result_index + 1
      
      cat(
        "results:",
        nrow(res)
      )
    }
    
    cat("\n")
  }
  
  #
  # Combine all results.
  #
  
  if (length(results) == 0) {
    return(NULL)
  }
  
  dplyr::bind_rows(results)
}
