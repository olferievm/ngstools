


filter_instrument_range <- function(ff,
                                    channels = c("FSC-A",
                                                 "SSC-A")) {
  
  p <- pData(parameters(ff))
  e <- exprs(ff)
  
  keep <- rep(TRUE, nrow(e))
  
  for (ch in channels) {
    
    i <- match(ch, p$name)
    
    keep <- keep &
      e[, ch] >= p$minRange[i] &
      e[, ch] <= p$maxRange[i]
  }
  
  ff[keep, ]
}