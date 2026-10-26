#' @title Account for fallback and re-ascension
#' 
#' @description
#'Applies the fallback and re-ascension rates to account for fallback in the abundance estimates. Differential rates are applied to adults and jacks
#'release above and below Lower Granite Dam as juveniles. The fallback rate is assumed to be zero for wild and unassigned hatchery-origin fish. 
#' 
#' @param indiv A dataframe with individual-level data. 
#' 
#' @param fallback A dataframe with the fallback rates for adults and jacks released above and below Lower Granite Dam. 
#' 
#' @param jack A column assigning fish as jacks or adults. Must be coded identically to the \code{fallback} df. 
#' 
#' @param origin A column specifying origin (i.e., hatchery, wild or unassigned hatchery). Hatchery-origin fish must be coded as \code{hatchery}. 
#' 
#' @param rel_lgr A column specifiying whether the fish was released above or below Lower Granite Dam as a juvenile. Must be coded identically to the \code{fallback} df. 
#' 
#' @param count_col The current count column, to be used for the expansion. 
#' 
#' @return a dataframe with an additional column \code{fallback_expand}.
#' 
#' @export


apply_fallback <- function(indiv, 
                           fallback, 
                           jack = "jack", 
                           origin = "origin", 
                           rel_lgr = "RelLGR",
                           count_col = "gen_expand") {
  
  join_var <- setNames(c(jack, rel_lgr), c(jack, rel_lgr))
  
  fallback_account <- indiv %>%
    dplyr::left_join(fallback, by = join_var) %>%
    dplyr::mutate(fallback_rate = as.numeric(fallback_rate)) %>%
    dplyr::mutate(fallback_rate = 
                    dplyr::if_else(.data[[origin]] != "hatchery", 0, fallback_rate)) %>%
    dplyr::mutate(fallback_expand = .data[[count_col]]*(1-fallback_rate))
  
  fallback_account
  
}
