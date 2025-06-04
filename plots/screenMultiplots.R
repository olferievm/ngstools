#' Generate multiple plots from a plotparameterized Excel sheet
#'
#' This script defines a function `screenMultiplot` that takes a data frame and a plotparameter
#' table (read from an Excel file), and produces plots (boxplots or scatterplots) accordingly.
#' Helper functions are included for parsing the plotparameter strings.

library(dplyr)
library(stringr)
library(ggpubr)

# Main plotting function
screenMultiplot <- function(x, plotparam, path) {
  # x: data.frame with values to plot
  # plotparam: data.frame with plot plotparameters (each row = one plot)
  # path: directory to save output images

#  x <- tmp0
#  plotparam <-  bxp
#  
  plotparam <- plotparam %>%
    dplyr::mutate(across(where(is.character), ~ gsub("&#10;", "\n", .))) %>%
    dplyr::mutate(
      xlab = ifelse(is.na(xlab), "", xlab),
      ylab = ifelse(is.na(ylab), "", ylab),
      comparisons = parse_str_to_list(comparisons),
      order = parse_order_to_list(order),
      fillcolors = parse_colors_to_str(fillcolors),
      label.x = suppressWarnings(as.numeric(label.x)),
      label.y = suppressWarnings(as.numeric(label.y))
    )
  
 
  
  for (i in seq_len(nrow(plotparam))) {

#    i <- 11
    g <- NULL
    
    df <- x %>%
      dplyr::filter(!is.na(.data[[plotparam$xpar[i]]]),
                    !is.na(.data[[plotparam$ypar[i]]]))
    
    if (plotparam$type[i] == 'boxplot') {
      
      g <- readyboxplot(df = df, plot_param = plotparam[i, ])
      
    } else if (plotparam$type[i] == 'scaterplot') {

      g <- readyscaterplot(df = df, plot_param = plotparam[i, ])
    }
    
    if (!is.null(g)) {
      ggsave(
        filename = file.path(path, sprintf("gg%s_%s_per_%s_%s.png",
                                           plotparam$type[i], plotparam$ypar[i], plotparam$xpar[i], plotparam$name[i])),
        plot = g,
        width = plotparam$fig.width[i],
        height = plotparam$fig.height[i]
      )
      cat(plotparam$type[i], " ", plotparam$name[i], "saved\n")
    }
  }
}

# Helper function to parse comparison string into list of lists
parse_str_to_list <- function(x) {
  lapply(x, function(mystring) {
    if (is.na(mystring) || mystring == "") return(NA)
    stringr::str_split_1(mystring, pattern = "\\n") %>%
      trimws() %>% stringr::str_split(., pattern = ",\\s*") %>% .[. !=""] %>% lapply(.,trimws)
  })
}

# Helper to parse order string into vector
parse_order_to_list <- function(x) {
  lapply(x, function(vct) {
    if (is.na(vct) || vct == "") return(NULL)
    stringr::str_split_1(vct, pattern = ",\\s*") %>% trimws()
  })
}



# Helper to parse color assignment string into named vector
parse_colors_to_str <- function(x) {
  lapply(x, function(colvect) {
    if (is.na(colvect) || colvect == "") return(NULL)
    colvect <- stringr::str_split(
      colvect, pattern = "\n")[[1]] %>%
      stringr::str_split(., pattern = "=")
    out <- sapply(colvect,'[[',2)
    names(out) <- sapply(colvect,'[[',1)
    return(out)
  })
}





# Boxplot generator
readyboxplot <- function(df, plot_param) {

  g <- ggboxplot(df,
                 x = plot_param$xpar[[1]],
                 y = plot_param$ypar[[1]],
                 bxp.errorbar = TRUE,
                 add = "jitter",
                 add.params = list(color = "black"),
                 xlab = plot_param$xlab[[1]],
                 ylab = plot_param$ylab[[1]],
                 fill = plot_param$fill[[1]],
                 palette = plot_param$fillcolors[[1]],
                 order = plot_param$order[[1]],
                 facet.by = plot_param$facet.by[[1]]
  )
  
  if (!any(is.na(plot_param$comparisons[[1]]))) {
    g <- g + stat_compare_means(comparisons = plot_param$comparisons[[1]])
  }
  return(g)
}

# Scatterplot generator
readyscaterplot <- function(df, plot_param) {
  # for test plot_param <- plotparam[i, ]
  g <- ggscatter(df,
                 x = plot_param$xpar[[1]],
                 y = plot_param$ypar[[1]],
                 add = "reg.line",
                 conf.int = TRUE,
                 xlab = plot_param$xlab[[1]],
                 ylab = plot_param$ylab[[1]],
                 color = plot_param$color[[1]],
                 palette = plot_param$fillcolors[[1]],
                 facet.by = plot_param$facet.by[[1]]
  )
  
  if (!any(is.na(plot_param$lm_stat[[1]]))) {
    g <- g + stat_cor(
      aes(color = plot_param$lm_stat[[1]]),
      label.x = plot_param$label.x[[1]],
      label.y = plot_param$label.y[[1]]
    )
  }
  return(g)
}
