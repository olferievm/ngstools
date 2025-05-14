#' Test for outliers using IQR method
#'
#' Identifies which values in a numeric vector are outliers using the IQR rule.
#'
#' @param x A numeric vector.
#' @param threshold A numeric multiplier for the IQR (default is 1.5, as in Tukey's rule).
#' @return A logical vector indicating which values are outliers (TRUE = outlier).
#' @examples
#' x <- c(1, 2, 3, 4, 5, 100)
#' test_if_outlier(x)
test_if_outlier <- function(x, threshold = 1.5) {
  if (!is.numeric(x)) stop("Input x must be numeric")
  
  lower <- quantile(x, 0.25, na.rm = TRUE) - threshold * IQR(x, na.rm = TRUE)
  upper <- quantile(x, 0.75, na.rm = TRUE) + threshold * IQR(x, na.rm = TRUE)
  
  return(x < lower | x > upper)
}
