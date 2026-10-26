#' @title Apply a trap-outage expansion
#' 
#' @description
#' Uses the window count at Lower Granite Dam and composition proportions from the five days preceding and following the trap outage to estimate abundance and 
#' composition during the trap outage period. Composition is determined by release site (i.e., final PBT group), sex, jack vs. adult, age, and origin (i.e., unassigned or hatchery-origin).
#' 
#' This function may need updating to ensure that ages of wild fish are assigned correctly.  
#' 
#' @param beg_outage The first date of the trap outage period, specified as \code{"mm/dd/yyyy"}.
#' 
#' @param end_outage The last date of the trap outage period, specified as \code{"mm/dd/yyyy"}.
#' 
#' @param window_count The window count during the trap outage period.
#' 
#' @param indiv The dataframe with individual-level data.
#' 
#' @param count_col The current count column, to be used for the expansion. 
#' 
#' @param agePBT The column specifying ages for fish assigned to PBT populations.
#' 
#' @param rel_site The column specifying the release site or final release group for fish assigned to a PBT population.
#' 
#' @param sex_col The column specifying the sex of the fish.
#' 
#' @param jack_col The column classifying fish as adults or jacks
#' 
#' @return A dataframe 
#' 
#' @export


trap_outage_expansion <- function(beg_outage, 
                                  end_outage, 
                                  window_count, 
                                  indiv, 
                                  count_col = "pbt_expand",
                                  agePBT = "agePBT", 
                                  rel_site = "RelSite", 
                                  sex_col = "Sex", 
                                  jack_col = "jack") {
  
  before_outage <- as.Date(beg_outage, format = "%m/%d/%Y") - 5
  after_outage <- as.Date(end_outage, format = "%m/%d/%Y") + 5
  
  #release group and brood year proportions estimates from period immediately before and after outage
  outage_count <- indiv %>%
    dplyr::filter(DateSampled >= before_outage & DateSampled <= after_outage) %>%
    dplyr::group_by(.data[[agePBT]], .data[[rel_site]]) %>%
    dplyr::summarise(count = sum(.data[[count_col]]), .groups = "drop_last") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(prop = count/sum(count)) %>%
    dplyr::mutate(outage_est = prop * window_count) %>%
    dplyr::select(!count)
  
  #sex proportions from entire run
  sex_ratio <- indiv %>%
    dplyr::mutate(Sex = 
                    dplyr::if_else(.data[[jack_col]] == "jack", "jack", .data[[sex_col]])) %>%
    dplyr::group_by(.data[[agePBT]], .data[[rel_site]], .data[[sex_col]]) %>%
    dplyr::summarise(count = sum(fallback_expand), .groups = "drop_last") %>%
    dplyr::mutate(sex_prop = count/sum(count)) %>%
    dplyr::select(!count)
  
  trap_outage_count <- outage_count %>%
    dplyr::left_join(sex_ratio, by = c(agePBT, rel_site)) %>%
    dplyr::mutate(outage_count = outage_est * sex_prop) %>%
    dplyr::select(all_of(c(agePBT, rel_site)), Sex, outage_count)
  
  release_group_locations <- indiv %>%
    dplyr::select(.data[[rel_site]], RelLGR) %>% distinct()
  
  trap_outage_count <- trap_outage_count %>%
    dplyr::mutate(origin = 
                    dplyr::if_else(.data[[rel_site]] == "Unassigned", "wild", "hatchery")) %>%
    dplyr::mutate(jack = 
                    dplyr::if_else(Sex == "jack", "jack", "adult")) %>%
    dplyr::left_join(release_group_locations, by = rel_site) %>%
    apply_fallback(fallback = fallback, 
                   count_col = "outage_count") %>%
    dplyr::select(all_of(c(agePBT, rel_site, sex_col)), fallback_expand) %>%
    dplyr::rename(outage_count = fallback_expand) 
  
}
