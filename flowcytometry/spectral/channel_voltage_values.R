#' Extract parameter voltages and metadata from a flow cytometry object
#'
#' Extracts selected parameter metadata from the FCS description of a
#' \code{flowFrame}. The returned information includes the parameter name,
#' detector voltage, channel label, and fluorochrome annotation.
#'
#' The function is primarily intended for extracting detector voltage
#' information for fluorescence or spectral channels from FCS metadata.
#'
#' @param ff A \code{flowFrame} containing FCS parameter metadata.
#'
#' @param channels Optional character vector specifying parameter names to
#'   retain. Parameter names are matched against the \code{$PnN} field in
#'   the FCS description. If \code{NULL}, metadata for all parameters are
#'   returned.
#'
#' @return A data frame with one row per FCS parameter and the following
#'   columns:
#' \describe{
#'   \item{\code{N}}{Parameter name from the \code{$PnN} FCS keyword.}
#'   \item{\code{V}}{Detector voltage from the \code{$PnV} keyword.}
#'   \item{\code{L}}{Parameter label from the \code{$PnL} keyword.}
#'   \item{\code{Fluor}}{Fluorochrome annotation from the \code{$PnF}
#'   keyword.}
#' }
#'
#' Missing FCS keywords are returned as \code{NA}.
#'
#' @details
#' FCS parameter metadata are stored in the \code{description} slot of a
#' \code{flowFrame}. Parameter-specific keywords follow the FCS naming
#' convention \code{$P1N}, \code{$P1V}, \code{$P1L}, \code{$P1F}, and so on.
#'
#' The \code{channels} argument is useful for extracting voltages only for
#' fluorescence or spectral parameters used in downstream analysis.
#'
#' @examples
#' \dontrun{
#' # Extract metadata for all parameters
#' channel_voltage_values(ff)
#'
#' # Extract metadata for selected fluorescence channels
#' channel_voltage_values(
#'   ff,
#'   channels = c("B488-A", "YG561-A", "R637-A")
#' )
#' }



channel_voltage_values <- function(ff, channels = NULL) {
  
  desc <- ff@description
  
  npar <- as.integer(desc[["$PAR"]])
  
  # Helper function returning NA for missing FCS keywords
  get_keyword <- function(i, keyword) {
    
    value <- desc[[paste0("$P", i, keyword)]]
    
    if (is.null(value)) {
      NA_character_
    } else {
      as.character(value)
    }
  }
  
  out <- data.frame(
    N = vapply(
      seq_len(npar),
      get_keyword,
      character(1),
      keyword = "N"
    ),
    V = vapply(
      seq_len(npar),
      get_keyword,
      character(1),
      keyword = "V"
    ),
    L = vapply(
      seq_len(npar),
      get_keyword,
      character(1),
      keyword = "L"
    ),
    Fluor = vapply(
      seq_len(npar),
      get_keyword,
      character(1),
      keyword = "F"
    ),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(channels)) {
    out <- out[out$N %in% channels, , drop = FALSE]
  }
  
  rownames(out) <- NULL
  
  out
}

channel_voltage_values(ff)
