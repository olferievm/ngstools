#' Arcsinh transformation for Flow CYtometry
#'
#' @param x Value to transform
#' @param a shift parameter
#' @param b scale parameter 5 for CyTOF (mass cytometry) data and 150 to 3000 (Cytek) for fluorescent or spectral flow cytometry data
#' @param c additional offset
#' @export
#' 
#' 

arcsinhTransform <- function(x,  a = 0, b = 150, c = 0) {
  
  x  <-  asinh(a + x / b) + c

}

