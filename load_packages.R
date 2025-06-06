#' Load and Install CRAN, Bioconductor, and GitHub Packages
#'
#' This function checks for and installs any missing packages from CRAN, Bioconductor,
#' or GitHub (or other remote repos via `remotes::install_*`), and then loads them.
#'
#' @param pkgs A character vector of all package names (CRAN and Bioconductor).
#' @param bioc_pkgs A character vector of package names known to come from Bioconductor.
#' @param github_pkgs A named character vector of GitHub repo strings (e.g., "user/repo"), 
#'                    where the names are the resulting package names.
#' @param quietly Logical; whether to suppress loading messages. Default is TRUE.
#'
#' @return Invisibly returns TRUE when all packages are loaded.
#' @examples
#' \dontrun{
#'   cran_pkgs <- c("dplyr", "ggplot2")
#'   bioc_pkgs <- c("GSVA")
#'   github_pkgs <- c("plotcor" = "yourgithub/plotcor")
#'   load_all_packages(pkgs = c(cran_pkgs, bioc_pkgs, names(github_pkgs)),
#'                     bioc_pkgs = bioc_pkgs,
#'                     github_pkgs = github_pkgs)
#' }
#' @export
load_packages <- function(pkgs,
                              bioc_pkgs = NULL,
                              github_pkgs = NULL,
                              quietly = TRUE) {
  # Install BiocManager if needed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # Install remotes if needed (for GitHub installs)
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  
  # Identify which pkgs are CRAN (excluding BioC and GitHub)
  cran_pkgs <- setdiff(pkgs, c(bioc_pkgs, names(github_pkgs)))
  
  # Load CRAN packages
  for (pkg in cran_pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = quietly)) {
      message("Installing CRAN package: ", pkg)
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
  
  # Load Bioconductor packages
  for (pkg in bioc_pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = quietly)) {
      message("Installing Bioconductor package: ", pkg)
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
      library(pkg, character.only = TRUE)
    }
  }
  
  # Load GitHub (or remote) packages
  for (pkg in names(github_pkgs)) {
    if (!require(pkg, character.only = TRUE, quietly = quietly)) {
      repo <- github_pkgs[[pkg]]
      message("Installing GitHub package: ", pkg, " from ", repo)
      remotes::install_github(repo, upgrade = "never", quiet = quietly)
      library(pkg, character.only = TRUE)
    }
  }
  
  invisible(TRUE)
}
