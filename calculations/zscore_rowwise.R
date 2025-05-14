# Row-wise z-score using scale()
zscore_rowwise <- function(x) {
  t(scale(t(x)))
}
