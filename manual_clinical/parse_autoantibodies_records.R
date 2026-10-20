parse_autoantibodies_records <- function(data, column = "autoantibodies") {
  
  # Note: The list of unique reported tests was collected and alphabetically sorted
  # The critical issue is that some names are prefixes of other names:
  # aCL → aCL_IgA, aCL_IgG, aCL_IgM
  # b2GP1 → b2GP1_IgA, etc.
  # DAT → DAT_COLD, DAT_WARM
  # Reticulin → Reticulin_IgA
  # SSA → SSA52, SSA60
  
  # Preserve the distinction between a general test and a specific test.
  
  antibody_names <- c(
    "aCL", "aCL_IgA", "aCL_IgG", "aCL_IgM",
    "ANA", "ANCA", "cANCA", "pANCA", "anti-C3_IgG", "anti-Gliadin", "anti-P",
    "aPL", "ASCA", "ASLO", "ATA", "ATG",
    "b2GP1", "b2GP1_IgA", "b2GP1_IgG", "b2GP1_IgM",
    "C1q", "CCP", "Centromere", "Chromatin", "CMV_IgG", "CMV_IgM", "Cryoglobulin",
    "DAT", "DAT_COLD", "DAT_WARM", "DRVVT", "dsDNA", "ENA", "EMA_IgA",
    "F-actin", "FTA", "Jo1", "GBM", "Gliadin_IgA", "Gliadin_IgG",
    "HepB_SA", "HepC", "Histone",
    "HLAB27", "hypergammaglobulinemia", "Ku",  "LAC", "Lyme", "MDA5",
    "MPO", "Mi2", "Myocardial_Ab", "Mitochondrial", "Nxp2",  "pANCA",
    "Phosphotodylserine_IgG", "Phosphotodylserine_IgM",
    "PL7", "PL12","PR", "PR3", "RBP",
    "Reticulin", "Reticulin_IgA", "RF", "RibosomalP",
    "RNAPol3", "RNP", "RPR", "SCL70", "Sm", "SmRNP", "SSA", "SSA52",
    "SSA60", "SSB", "Srp", "OJ", "EJ", "TG", "TPO", "TTG", "Tif1y", "VDRL"
  )
  
 # "   Ku-; TG-; TPO-; Mi2-; PL12-; PL7-; Srp-; Jo1-; OJ-; EJ-; b2GP1_IgG-; b2GP1_IgM-; b2GP1_IgA-;  aCL_IgG-;   aCL_IgM-;  aCL_IgA-;";
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
  #column <- 'Serology'
  #data <- serol_df
  #rm(data)
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