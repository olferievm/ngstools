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
    set4Group = c(
      green ="#4DAF4A",
      blue = "#377EB8",
      orange = "#FF7F00",
      red = "#E41A1C"
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
      blue1 = "#0803A1",
      blue2 = "#1514A6",
      blue3 = "#225EA8",
      blue4 = "#74A9CF",
      gray1 = "#BDC9E1",
      gray2 = "#F1EEF6",
      yellow1 = "#FDCC8A",
      orange1 = "#FC8D59",
      orange2 = "#D7301F",
      red1 = "#C11C38",
      red2 = "#b01919"
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
    ),
    descr_BlOr = c(
      blue = "#377EB8",
      orange = "#FF7F00"),
    descr_NavyRed = c(
      navy = "navy",
      darkred = "darkred"),
    gr2ColBlRd = c(
      navy = '#172767',
      red = "#b80422"),
    gr3ColBlYlRd = c(
      navy = '#172767',
      orange = "#ee9b43",
      red = "#b80422"),
    grExterCol = c(
      yellow1 = "#ffec9d",
      yellow2 = "#fac881",
      yellow3 = "#f4a464",
      orange1 = "#e87444",
      orange2 = "#d9402a",
      red1 = "#bf2729",
      red2 = "#912534",
      red3 = "#64243e")
  )
  return(palettes)
}