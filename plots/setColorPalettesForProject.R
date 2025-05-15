#' Set Color Palettes for Project
#'
#' This function defines a list of custom color palettes commonly used in visualizations.
#' It includes a modified Set1 palette, gradient blues, a yellow-pink-violet trio,
#' a blue-white-yellow gradient, and a palette inspired by the Journal of Clinical Oncology.
#'
#' @return A named list containing the following color palettes:
#' \describe{
#'   \item{setCol}{Modified Set1 palette from RColorBrewer with colors like green, blue, red, etc.}
#'   \item{grPal}{Gradient palette of blues, grays, and warm tones.}
#'   \item{YlPkVi}{Yellow-pink-violet palette.}
#'   \item{myBlWtYl}{Blue-white-yellow gradient palette.}
#'   \item{jco}{Journal of Clinical Oncology-inspired palette.}
#' }
#' @examples
#' palettes <- setColorPalettesForProject()
#' palettes$setCol["green"]
#' 
#' @export

setColorPaletteForProject <- function() {
  palettes <- list(
    setCol = c(
      green  = "#4DAF4A",
      blue   = "#377EB8",
      orange = "#FF7F00",
      red    = "#E41A1C",
      purple = "#984EA3",
      yellow = "#FFFF33",
      brown  = "#A65628",
      pink   = "#F781BF"
    ),
    grPal = c(
      deepblue  = "#225EA8",
      skyblue   = "#74A9CF",
      paleblue  = "#BDC9E1",
      lightgray = "#FAFAFA",
      peach     = "#FDCC8A",
      coral     = "#FC8D59",
      brickred  = "#D7301F"
    ),
    YlPkVi = c(
      orange = "#FF7F00",
      pink   = "#FF40FF",
      blue   = "#011993"
    ),
    myBlWtYl = c(
      "#0803A1", "#1514A6", "#225EA8", "#74A9CF", "#BDC9E1",
      "#F1EEF6", "#FDCC8A", "#FC8D59", "#D7301F", "#C11C38", "#b01919"
    ),
    jco = c(
      blue       = "#0073C2",
      yellow     = "#EFC000",
      gray       = "#868686",
      red        = "#CD534C",
      lightblue  = "#7AA6DC",
      darkblue   = "#003C67",
      brown      = "#8F7700",
      darkgray   = "#3B3B3B",
      darkred    = "#A73030",
      slateblue  = "#4A6990"
    )
  )
  return(palettes)
}