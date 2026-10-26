#' @title Assign ages to unassigned hatchery-origin fish
#' 
#' @description
#' Uses an age-length key to assign ages to unassigned hatchery-origin fish. 
#' 
#' @param hat_unk A dataframe with the total estimated number of unassigned hatchery-origin fish, grouped by sex, release group, length bin, and jack vs adult. 
#' 
#' @param age_length_key An age-length key for natural origin fish for the specified sex.
#' 
#' @param sex A character specifying which sex to assign ages to. 
#' 
#' @param pool_FL A column specifying the length bin for each fish. The length bins in the dataframe must match the length bins in the age-length key.
#' 
#' @param count The current count column.  
#' 
#' @return A dataframe 
#' 
#' @export

assign_hat_unk_age <- function(hat_unk, 
                               age_length_key, 
                               sex = "F",
                               pool_FL = "pool_FL", 
                               count = "count") {
  
  ages <- colnames(age_length_key)[2:dim(age_length_key)[2]]
  
  hat_unk %>%
    dplyr::filter(Sex == sex) %>%
    dplyr::left_join(age_length_key, by = pool_FL) %>%
    tidyr::pivot_longer(cols = 
                   dplyr::all_of(ages), names_to = "age", values_to = "prop") %>%
    dplyr::mutate(count_by_age = prop * .data[[count]]) %>%
    dplyr::mutate(rearing = "subyearling") %>%
    dplyr::select(-count, -prop, -agePBT)  %>%
    dplyr::rename(count = count_by_age) %>%
    dplyr::mutate(age = as.numeric(age)) %>%
    dplyr::mutate(brood_year = run_yr - age)
}