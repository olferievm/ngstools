# Those are small functions used for data processing.
#'@param keep_one - function will return the first non-NA element, otherwise NA 
#'@param catch_wrong - Replace Excel string "N/A" with NA
#'@param clean_term -Remove parenthetical notes and trim whitespace
#'@param split_term Remove parenthetical notes, trim white space, split text by pattern
#'@param count_diff Return x - y only if both are non-NA scalars
#'@param
#'@param

# Return the first non-NA element, otherwise NA
keep_one <- function(x) {

  x <- na.omit(x)
  if (length(x)) x[[1]] else NA
}

# Summary table:
# Calculate cases with AVN and AVE other adverse events based on all visit.
#keep_one <- function(x){ if(all(is.na(x))){ return(NA)}
#  unique(x[!is.na(x)])[1]}

# Replace Excel string "N/A" with NA
catch_wrong <- function(x) {
  replace(x, x == "N/A", NA)
}


# Remove parenthetical notes and trim whitespace
clean_term <- function(x) {
  x <- gsub("\\s*\\([^)]*\\)", "", x)
  trimws(x)
}

# Remove parenthetical notes, trim white space, split text by pattern
split_term <- function(x, sep, fixed = TRUE) {
#'@param x - character vector
#'@param sep - pattern to split "\n", ";"
#'@param fixed - use character or 'regex' pattern 
#  split_term(x, "\n")
#  split_term(x, ";")
  strsplit(clean_term(x), sep, fixed = fixed)
}

split_term_newline <- function(x) {
  strsplit(clean_term(x), "\n", fixed = TRUE)
}


# Split cleaned text by semicolon
split_term_semicolon <- function(x) {
  strsplit(clean_term(x), ";", fixed = TRUE)
}


# Return x - y only if both are non-NA scalars
count_diff <- function(x, y) {
  if (length(x) != 1L || length(y) != 1L || is.na(x) || is.na(y)) {
    NA_real_
  } else {
    x - y
  }
}


