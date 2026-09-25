glmmLasso_variable_selection <- function(
    data,
    var_names,
    outcome,
    covariates = NULL,
    random_effect = "ID",
    lambda_seq = 10^seq(-1, 3, length.out = 30),
    nfolds = 10,
    seed = 37912,
    final.re = TRUE
) {
  
  set.seed(seed)
  
  # ------------------------------------------------------------
  # Basic checks
  # ------------------------------------------------------------
  stopifnot(outcome %in% names(data))
  stopifnot(all(var_names %in% names(data)))
  
  if (!is.null(covariates)) {
    stopifnot(all(covariates %in% names(data)))
  }
  
  stopifnot(random_effect %in% names(data))
  
  if(class(data[[random_effect]]) != "factor"){
     data[[random_effect]] <- factor(data[[random_effect]])
  }
  
  # Remove outcome from candidate predictors
  var_names <- setdiff(var_names, outcome)
  
  # ------------------------------------------------------------
  # Fixed-effect formula
  # ------------------------------------------------------------
  rhs <- c(covariates, var_names)
  
  fix_formula <- as.formula(
    paste(outcome, "~", paste(rhs, collapse = " + "))
  )
  
  # ------------------------------------------------------------
  # Subject-level folds
  # ------------------------------------------------------------
  IDs <- unique(data[[random_effect]])
  
  if (nfolds > length(IDs)) {
    nfolds <- length(IDs)
  }
  
  fold_id <- sample(rep(seq_len(nfolds), length.out = length(IDs)))
  names(fold_id) <- IDs
  
  # Index for glmmLasso:
  # positions of candidate variables in the data
  index <- which(colnames(data) %in% var_names)
  
  # ------------------------------------------------------------
  # Random-effect specification
  # ------------------------------------------------------------
  rnd_formula <- list()
  rnd_formula[[random_effect]] <- ~1
  
  # ------------------------------------------------------------
  # Cross-validation
  # ------------------------------------------------------------
  cv_results <- vector("list", length(lambda_seq))
  coef_results <- vector("list", length(lambda_seq))
  
  for (i in seq_along(lambda_seq)) {
    #i <- 1
    lambda <- lambda_seq[i]
    
    fold_mse <- numeric(nfolds)
    
    coef_ls <- vector("list", nfolds)
    
    for (k in seq_len(nfolds)) {
      #k = 2
      cat("lambda =",lambda,"; k = ",k,"\n")
      test_ids <- IDs[fold_id == k]
      
      train <- data[!data[[random_effect]] %in% test_ids, , drop = FALSE]
      test  <- data[ data[[random_effect]] %in% test_ids, , drop = FALSE]
      
      #lambda = 1
      fit <- tryCatch(
        glmmLasso::glmmLasso(
          fix = fix_formula,
          rnd = rnd_formula,
          data = train,
          lambda = lambda,
          final.re = final.re,
          control = list(
            index = index,
            standardize = TRUE,
            center = TRUE
          )
        ),
        
        error = function(e) NULL
      )
      
      if (is.null(fit)) {
        fold_mse[k] <- NA_real_
        next
      }
      
      # Fixed-effect prediction.
      # This is appropriate because test subjects are unseen IDs.
      X_test <- model.matrix(fix_formula, data = test)
      
      beta <- fit$coefficients
      
      coef_ls[[k]] <- beta %>% as.matrix(.) %>% as.data.frame(.) %>% tibble::rownames_to_column(., var = 'feature')
      coef_ls[[k]]$lambda <- lambda
      coef_ls[[k]]$n_fold <- k
      
      # Match model-matrix columns to fitted coefficients
      common <- intersect(colnames(X_test), names(beta))
      
      if (length(common) == 0) {
        fold_mse[k] <- NA_real_
        next
      }
      
      pred <- as.numeric(
        X_test[, common, drop = FALSE] %*% beta[common]
      )
      
      fold_mse[k] <- mean(
        (test[[outcome]] - pred)^2,
        na.rm = TRUE
      )
    }
    
    coef_ls <- coef_ls[!sapply(coef_ls,is.null)]
    coef_results[[i]] <- do.call(rbind, coef_ls)
    
    cv_results[[i]] <- data.frame(
      lambda = lambda,
      mse = mean(fold_mse, na.rm = TRUE),
      sd_mse = sd(fold_mse, na.rm = TRUE),
      n_valid_folds = sum(is.finite(fold_mse))
    )
  }
  
  cv_results <- do.call(rbind, cv_results)
  coef_results <- do.call(rbind, coef_results)
  # ------------------------------------------------------------
  # Select lambda
  # ------------------------------------------------------------
  best_row <- which.min(cv_results$mse)
  
  lambda_min <- cv_results$lambda[best_row]
  
  # ------------------------------------------------------------
  # Fit final model using all subjects
  # ------------------------------------------------------------
  final_fit <- glmmLasso::glmmLasso(
    fix = fix_formula,
    rnd = rnd_formula,
    data = data,
    lambda = lambda_min,
    final.re = final.re,
    control = list(
      index = index,
      standardize = TRUE,
      center = TRUE
    )
  )
  
  # ------------------------------------------------------------
  # Selected variables
  # ------------------------------------------------------------
  coefficients <- final_fit$coefficients
  
  selected <- var_names[
    var_names %in% names(coefficients) &
      coefficients[var_names] != 0
  ]
  
  list(
    model = final_fit,
    cv_results = cv_results,
    lambda.min = lambda_min,
    selected = selected,
    coefficients = coefficients,
    cv = cv_results,
    folds = fold_id,
    formula = fix_formula,
    index = index,
    coef_results = coef_results
  )
}