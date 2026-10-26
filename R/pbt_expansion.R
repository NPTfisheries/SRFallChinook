#' @title Apply an expansion for the PBT tagging rate. 
#' 
#' @description
#' Applies the expansion for the the PBT tagging rate by population. Each hatchery-assigned group is expanded based on the corresponding PBT tagging rate. This expansion adds to the 
#' number of hatchery-assigned fish, and therefore requires that unassigned fish are substracted. The total number of ad-intact and ad-clipped fish added to the hatchery-assigned 
#' group is determined, and the corresponding number of fish are removed from the ad-intact and ad-clipped unassigned groups. 
#' 
#' NOTE: To match the results of the 2024 and 2025 run reconstruction analyses, the expansion is applied as an additive expansion per row. However,
#' this is likely appropriate only if the PBT tag rate expansion were the first expansion applied. Otherwise, the expansion should be applied as a
#' multiplicative expansion.
#' 
#' @param indiv A dataframe with individual-level data. 
#' 
#' @param PbtTagRate A column with the PBT tag rates
#' 
#' @param PopPBT A column with the PBT population assignments. 
#' 
#' @param count_col The current count column, to be used for the expansion. 
#' 
#' @param Clip A column indicating if the fish was ad-clipped or ad-intact. 
#' 
#' @details A fish that was not assigned to a PBT population should be coded as \code{Unassigned} in the PopPBT column. 
#' 
#' @return a dataframe with an additional column \code{pbt_expand}
#' 
#' @export

pbt_expansion <- function(indiv, 
                          PbtTagRate = "PbtTagRate", 
                          count_col = "ng_expand", 
                          PopPBT = "PopPBT", 
                          Clip = "Clip") {
  
  pbt_expand_tbl <- indiv %>%
    #set PBT tag rate for unassigned fish to 1
    dplyr::mutate(PbtTagRate = 
                    dplyr::if_else(.data[[PopPBT]] == "Unassigned", 1, .data[[PbtTagRate]])) %>%
    #divide count by tag rate
    dplyr::mutate(pbt_expand = .data[[count_col]]/PbtTagRate) %>%
    dplyr::group_by(.data[[Clip]], .data[[PopPBT]]) %>%
    dplyr::summarise(pbt_expand = sum(pbt_expand), 
              previous_expand = sum(.data[[count_col]]), 
              diff = previous_expand - pbt_expand, .groups = "drop_last") %>%
    dplyr::group_by(Clip) %>%
    dplyr::mutate(exp_sum = sum(diff, na.rm = T), 
           diff = 
             dplyr::if_else(PopPBT == "Unassigned", exp_sum * -1, diff)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(prop_exp = diff/previous_expand)
  
  pbt_expand <- indiv %>%
    dplyr::left_join(pbt_expand_tbl %>%
                       dplyr::select(Clip, PopPBT, prop_exp), by = c("PopPBT", "Clip")) %>%
    dplyr::mutate(pbt_expand = .data[[count_col]] - prop_exp) %>%
    #mutate(pbt_expand = .data[[count_col]] * (1 - prop_exp)) %>%
    dplyr::select(-prop_exp)
  
  pbt_expand
}
