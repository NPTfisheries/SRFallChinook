#' @title Apply the genotype sampling rate expansion
#' 
#' @description
#' Applies the expansion for the genotype sampling rate at Lower Granite Dam by dividing the count of fish by the genotyping rate. 
#' 
#' @param count_col The current count column, to be used for the expansion. 
#' 
#' @param gen_rate A column containing the corresponding genotype sampling rate. 
#' 
#' @return a dataframe with an additional column \code{gen_rate_expand}
#' 
#' @export


gen_rate_expansion <- function(indiv, 
                               count_col = "pbt_expand", 
                               gen_rate = "GenRate") {
  gen_rate_expand <- indiv %>%
    dplyr::mutate(gen_expand = .data[[count_col]]/.data[[gen_rate]])
  
  gen_rate_expand
}
