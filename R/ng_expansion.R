#' @title Apply an expansion for samples that failed to genotype
#' 
#' @description
#' Applies the expansion for the samples that failed to genotype successfully. The expansion is differentially by groups, specified by \code{group_var}. The current run
#' reconstruction method applies this expansion by pooled fork length, sex, trap rate period, and by release site. Samples that did not genotype are then removed from the returned
#' dataframe so they are not included in future counts.   
#' 
#' @param indiv A dataframe with individual-level data. 
#' 
#' @param group_var A vector with the names of the grouping variables for the expansion
#' 
#' @param ng_col The column with the information on whether the sample genotyped successfully. In the current datasets, 
#' this is in the column specifying the release group the sample was assigned to. 
#' 
#' @param count_col The current count column, to be used for the expansion. 
#' 
#' @details A sample that failed to genotype must be coded as \code{NG}.  
#' 
#' @return a dataframe with an additional column \code{ng_expand}
#' 
#' @export

ng_expansion <- function(indiv, 
                         group_var = c("pool_FL", "Sex", "TrapRate", "RelSite"), 
                         ng_col = "RelSite",
                         count_col = lgr_expand) {
  
  pool_var <- dplyr::setdiff(group_var, ng_col)
  
  #calculate the rate of failed genotyping by length, sex, and trap rate
  ng_prop_tbl <- indiv %>%
    dplyr::group_by(across(
      dplyr::all_of(group_var))) %>%
    dplyr::summarise(count = sum(.data[[count_col]]), .groups = "drop_last") %>%
    dplyr::mutate(gen_prop = 1 - (sum(count[.data[[ng_col]] == "NG"])/sum(count))) %>%
    dplyr::ungroup() %>%
    dplyr::select(all_of(pool_var), gen_prop) %>%
    dplyr::distinct()
  
  #expand by the genotyping rate, proportional based on length, sex, and trap rate
  #set NG fish count to 0
  #remove NG fish from the indiv dataset
  ng_expand <- indiv %>%
    dplyr::full_join(ng_prop_tbl, by = pool_var) %>%
    dplyr::mutate(ng_expand = .data[[count_col]]/gen_prop) %>%
    dplyr::mutate(ng_expand = if_else(.data[[ng_col]] == "NG", 0, ng_expand)) %>%
    dplyr::filter(.data[[ng_col]] != "NG")
  
  ng_prop_tbl
  ng_expand
  
}  
