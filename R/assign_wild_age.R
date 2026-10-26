#' @title Assign ages to natural-origin fish
#' 
#' @description
#' Uses an age-length key to assign ages to wild fish. European ages are used in the age-length key are converted to a simple age, defined as the run year - the brood year. 
#' 
#' @param wild A dataframe with the total estimated number of natural-origin fish, grouped by sex, length bin, and jack vs adult. 
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

assign_wild_age <- function(wild,
                            age_length_key, 
                            sex = "F",
                            pool_FL = "pool_FL", 
                            count = "count") {
  scale_ages <- colnames(age_length_key)[2:dim(age_length_key)[2]]
  
  wild %>%
    dplyr::filter(Sex == sex) %>%
    dplyr::left_join(age_length_key, by = "pool_FL") %>%
    tidyr::pivot_longer(cols = all_of(scale_ages), names_to = "scale_age", values_to = "prop") %>%
    dplyr::mutate(count_by_age = prop * .data[[count]]) %>%
    dplyr::mutate(scale_age = 
                    dplyr::if_else(scale_age == "1", "1.0", scale_age)) %>%
    dplyr::mutate(age = as.numeric(substr(scale_age, 1, 1)) + as.numeric(substr(scale_age, 3, 3)) + 1) %>%
    dplyr::mutate(rearing = 
                    dplyr::if_else(substr(scale_age, 1, 1) == 0, "subyearling", "yearling")) %>%
    dplyr::select(-count, -prop, -scale_age, -agePBT) %>%
    dplyr::rename(count = count_by_age) %>%
    dplyr::mutate(brood_year = run_yr - age)
}
