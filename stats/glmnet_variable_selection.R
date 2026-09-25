glmnet_variable_selection <- function(
    data,
    var_names,
    outcome = NULL,
    covariates = NULL,
    alpha = 1,
    nfolds = 10,
    seed = 37912
) {
  
  # ----------------------------------------------------------
  # Check variables
  # ----------------------------------------------------------
  
  var_names <- intersect(var_names, names(data))
  
  if (length(var_names) == 0) {
    stop("None of the variables in 'var_names' were found in 'data'.")
  }
  
  if (!is.null(outcome) && !outcome %in% names(data)) {
    stop("Outcome variable not found in 'data'.")
  }
  
  if (!is.null(covariates)) {
    covariates <- intersect(covariates, names(data))
    
    if (length(covariates) == 0) {
      stop("None of the covariates were found in 'data'.")
    }
  }
  
  # Outcome cannot also be a predictor
  if (!is.null(outcome)) {
    var_names <- setdiff(var_names, outcome)
  }
  
  # ----------------------------------------------------------
  # Select complete cases
  # ----------------------------------------------------------
  
  vars_needed <- unique(
    c(outcome, var_names, covariates)
  )
  
  dat <- data %>%
    dplyr::select(dplyr::all_of(vars_needed)) %>%
    dplyr::filter(stats::complete.cases(.))
  
  if (nrow(dat) < nfolds) {
    stop("Number of complete cases is smaller than nfolds.")
  }
  
  # ----------------------------------------------------------
  # Outcome
  # ----------------------------------------------------------
  
  if (is.null(outcome)) {
    stop("'outcome' must be specified for glmnet variable selection.")
  }
  
  y <- dat[[outcome]]
  
  if (!is.numeric(y)) {
    stop("Outcome must be numeric for family = 'gaussian'.")
  }
  
  # ----------------------------------------------------------
  # Candidate predictors
  # ----------------------------------------------------------
  
  X_var <- as.matrix(
    dat[, var_names, drop = FALSE]
  )
  
  storage.mode(X_var) <- "numeric"
  
  # ----------------------------------------------------------
  # Covariates
  #
  # model.matrix() automatically converts factors into
  # dummy variables.
  # ----------------------------------------------------------
  
  if (!is.null(covariates)) {
    
    for (v in covariates) {
      dat[[v]] <- factor(dat[[v]])
    }
    
    covar_formula <- as.formula(
      paste(
        "~",
        paste(covariates, collapse = " + "),
        "-1"
      )
    )
    
    X_cov <- model.matrix(
      covar_formula,
      data = dat
    )
    
  } else {
    
    X_cov <- matrix(
      nrow = nrow(dat),
      ncol = 0
    )
  }
  
  # ----------------------------------------------------------
  # Combine covariates + variables
  # ----------------------------------------------------------
  
  X <- cbind(
    X_cov,
    X_var
  )
  
  # ----------------------------------------------------------
  # Penalty factors
  #
  # Covariates: 0 = never penalized
  # Candidate variables: 1 = penalized
  # ----------------------------------------------------------
  
  penalty.factor <- c(
    rep(0, ncol(X_cov)),
    rep(1, ncol(X_var))
  )
  
  # ----------------------------------------------------------
  # Cross-validated glmnet
  # ----------------------------------------------------------
  
  set.seed(seed)
  
  cvfit <- glmnet::cv.glmnet(
    x = X,
    y = y,
    family = "gaussian",
    alpha = alpha,
    nfolds = nfolds,
    standardize = TRUE,
    penalty.factor = penalty.factor,
    type.measure = "mse"
  )
  
  # ----------------------------------------------------------
  # Extract coefficients
  # ----------------------------------------------------------
  
  extract_coefficients <- function(s) {
    
    coef_mat <- as.matrix(
      stats::coef(cvfit, s = s)
    )
    
    coef_mat %>%
      as.data.frame() %>%
      tibble::rownames_to_column("variable") %>%
      dplyr::rename(coefficient = 2) %>%
      dplyr::filter(
        variable != "(Intercept)",
        coefficient != 0
      ) %>%
      dplyr::mutate(
        abs_coefficient = abs(coefficient),
        lambda = s
      ) %>%
      dplyr::arrange(dplyr::desc(abs_coefficient))
  }
  
  coef_min <- extract_coefficients("lambda.min")
  coef_1se <- extract_coefficients("lambda.1se")
  
  # ----------------------------------------------------------
  # Identify selected candidate variables
  # ----------------------------------------------------------
  
  selected_min <- coef_min %>%
    dplyr::filter(
      variable %in% var_names
    )
  
  selected_1se <- coef_1se %>%
    dplyr::filter(
      variable %in% var_names
    )
  
  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------
  
  list(
    data = dat,
    
    # Design matrix
    X = X,
    y = y,
    
    # glmnet object
    cvfit = cvfit,
    
    # Lambda values
    lambda_min = cvfit$lambda.min,
    lambda_1se = cvfit$lambda.1se,
    
    # All non-zero coefficients
    coefficients_lambda_min = coef_min,
    coefficients_lambda_1se = coef_1se,
    
    # Selected candidate variables only
    selected_lambda_min = selected_min,
    selected_lambda_1se = selected_1se,
    
    # Variables used
    var_names = var_names,
    covariates = covariates,
    
    # Scaling information
    means = cvfit$glmnet.fit$xm,
    scales = cvfit$glmnet.fit$xs
  )
}