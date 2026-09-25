parse_autoantibodies_records <- function(data, column = autoantibodies) {
  
  # Note: The list of unique reported tests was collected and alphabetically sorted
  # The critical issue is that some names are prefixes of other names:
  # aCL → aCL_IgA, aCL_IgG, aCL_IgM
  # b2GP1 → b2GP1_IgA, etc.
  # DAT → DAT_COLD, DAT_WARM
  # Reticulin → Reticulin_IgA
  # SSA → SSA52, SSA60
  
  # Preserve the distinction between a general test and a specific test.
  # b2GP1+ should produce 
  # b2GP1 = "+"
  # b2GP1_IgG = NA
  # b2GP1_IgM = NA
  # b2GP1_IgA = NA
  
  antibody_names <- c(
    "aCL", "aCL_IgA", "aCL_IgG", "aCL_IgM",
    "ANA", "ANCA", "anti-C3_IgG", "anti-Gliadin", "anti-P",
    "aPL", "ASCA", "ASLO", "ATA", "ATG",
    "b2GP1", "b2GP1_IgA", "b2GP1_IgG", "b2GP1_IgM",
    "C1q", "CCP", "Centromere", "Chromatin", "Cryoglobulin",
    "DAT", "DAT_COLD", "DAT_WARM", "DRVVT", "dsDNA", "ENA",
    "F-actin", "FTA", "GBM", "HepB_SA", "HepC", "Histone",
    "HLAB27", "hypergammaglobulinemia", "Jo1", "LAC", "Lyme",
    "MPO", "Myocardial_Ab", "pANCA", "PR", "PR3", "RBP",
    "Reticulin", "Reticulin_IgA", "RF", "RibosomalP",
    "RNAPol3", "RNP", "RPR", "SCL70", "Sm", "SSA", "SSA52",
    "SSA60", "SSB", "TG", "TPO", "TTG", "VDRL"
  )
  
  # Longer names first to prevent prefix matching
  antibody_names <- antibody_names[
    order(nchar(antibody_names), decreasing = TRUE)
  ]
  
  antibody_pattern <- paste(
    str_replace_all(
      antibody_names,
      "([.|()\\^{}+$*?\\[\\]\\\\])",
      "\\\\\\1"
    ),
    collapse = "|"
  )
  
  pattern <- paste0(
    "^(",
    antibody_pattern,
    ")(\\+/-|\\+|-|low|high)?$"
  )
  #data <- serol_df
  #column <- 'Serology'
  #z[4,]$.autoantibody_item
  parsed <- data %>%
    mutate(
      .autoantibody_item = str_split({{ column }}, ";\\s*")
    ) %>%
    unnest_longer(.autoantibody_item) %>%
    filter(
       !is.na(.autoantibody_item),
       .autoantibody_item != ""
     ) %>%
     mutate(
       .match = str_match(.autoantibody_item, pattern),
       .antibody = .match[, 2],
       .value = .match[, 3]
     ) %>%
     select(-.match, -.autoantibody_item)
  
  duplicates <- parsed %>%
    filter(!is.na(.antibody)) %>%
    group_by(ID, .antibody) %>%
    filter(n() > 1) %>%
    summarise(
      duplicated_items = paste(.antibody, collapse = "; "),
      .groups = "drop"
    )
  
  if (nrow(duplicates) > 0) {
    
    message(
      "\nDuplicated autoantibody items detected.\n",
      "The script will stop. Review the following records:\n"
    )
    
    print(duplicates, n = Inf)
    
    stop(
      "\nDuplicated autoantibody items detected. ",
      "Correct the source data before continuing."
    )
  }
  
  parsed %>%
     tidyr::pivot_wider(
       id_cols = c(ID, Serology),
       names_from = .antibody,
       values_from = .value,
       values_fill = NA_character_
     )
}