#' @title Create an age-length key
#' 
#' @description
#' Uses raw age data to create an age-length key. The age-length key is developed using pre-defined length bins. 
#' 
#' @param age A dataframe with ages for individual fish.
#' 
#' @param sex A character, defining whether the age-length key should be created for male or female fish.
#' 
#' @param sex_col A column specifying the sex of fish in the dataset. Must match the character defined in \code{sex}.
#' 
#' @param age_col A column with the age data
#' 
#' @param pool_FL A column specifying the length bin for each fish.  
#' 
#' @return A dataframe 
#' 
#' @export


create_age_length_key <- function(age, 
                                  sex = "M", 
                                  sex_col = "PhenotypicSex",
                                  age_col = "Age", 
                                  pool_FL = "pool_FL") {
  
  age_length_key <- age %>%
    dplyr::filter(.data[[sex_col]] == sex) %>%
    dplyr::group_by(.data[[age_col]], .data[[pool_FL]]) %>%
    count() %>% 
    dplyr::ungroup() %>%
    tidyr::complete(.data[[age_col]], .data[[pool_FL]], fill = list(n = 0)) %>%
    dplyr::group_by(.data[[pool_FL]]) %>%
    dplyr::mutate(prop = n/sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::select(dplyr::all_of(c(age_col, pool_FL)), prop) %>%
    tidyr::pivot_wider(names_from = dplyr::all_of(age_col), values_from = prop) 
  
  age_length_key
}
