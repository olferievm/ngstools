# Assume the data.frame includes columns
# PC01, PC02, PC03 ... etc numeric (or WGCNA vectors etc)
# RACE, GENDER, HGB, HCT ... etc factor or numeric
# Then find linear associations using the function below
# x = name of parameter ("HGB")
# pc = name of numerical parameter ("PC01", but can be any else)
# return statistics of linear model

#' @title PCA Covariate Correlation
#' @description
#'   Fits a linear model between a given variable and a principal component (PC) or other numerical parameter,
#'   extracts both the coefficient statistics and model fit statistics.
#' Assume the data.frame includes columns
#' PC01, PC02, PC03 ... etc numeric
#' RACE, GENDER, HGB, HCT ... etc factor or numeric
#' Then find linear associations using the function below
#' x = name of parameter ("HGB")
#' pc = name of numerical parameter ("PC01", but can be any else)
#' return statistics of linear model
#' @param x Character. Name of the independent variable (e.g., a gene or metadata column).
#' @param pc Character. Name of the principal component (e.g., "PC1", "PC2") to correlate against.
#' @param data data.frame with a choosen columns
#' @return A tibble with statistics of the linear model, including:
#'   \itemize{
#'     \item \code{PC}: The PC name.
#'     \item \code{term}: The covariate term (same as input \code{x}).
#'     \item \code{estimate}: Coefficient estimate.
#'     \item \code{std.error}: Standard error of the estimate.
#'     \item \code{stats}: t-statistic for the coefficient.
#'     \item \code{pval}: p-value for the coefficient.
#'     \item \code{r.squared}, \code{adj.r.squared}, \code{fit.statistic}, \code{fit.pvalue}, etc.
#'   }
#'
#' @details
#' This function assumes the existence of a data frame named `tmp` in the global environment.
#' You may consider passing the data frame as an argument for better portability.
#'
#' @importFrom broom tidy glance
#' @importFrom dplyr filter rename mutate select everything bind_cols
#' @examples
#' # Assuming 'tmp' contains columns "Gene1" and "PC1":
#' find_linear_associations("Gene1", "PC1")
#'


find_linear_associations <- function(x, pc, data) {
  # Construct a formula like PC1 ~ gene_expression
  formula <- as.formula(paste(pc, "~", x))
  
  # Fit linear model using the 'data'
  lm_fit <- lm(formula, data = data)
  
  # Extract and clean the model coefficient (excluding intercept), rename statistics
  coef_summary <- broom::tidy(lm_fit) %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::rename(
      stats = statistic,
      pval = p.value
    )
  
  # Extract model fit statistics and rename fields
  fit_summary <- broom::glance(lm_fit) %>%
    dplyr::rename(
      fit.statistic = statistic,
      fit.pvalue = p.value
    )
  
  # Combine coefficient stats and fit stats; add PC name and rearrange columns
  result <- dplyr::bind_cols(coef_summary, fit_summary) %>%
    dplyr::mutate(PC = pc) %>%
    dplyr::select(PC, dplyr::everything())
  
  return(result)
}
