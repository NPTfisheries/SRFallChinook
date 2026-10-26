#' @title Assigns fish origin
#' 
#' @description
#' Fish that were not assigned to a PBT population and are ad-intact are classified as \code{wild}; fish that were not assigned to a PBT population but 
#' are ad-clipped are classified as \code{hatchery_unknown}; fish that were assigned to a PBT population are classified as \code{hatchery}. 
#' 
#' @param Clip A column specifying whether a fish was ad-clipped or ad-intact. 
#' 
#' @param PopPBT A column specifying if a fish was assigned to a PBT population. 
#' 
#' @details Ad-clipped fish must be specified as \code{AD} and ad-intact fish must be specified as \code{AI}. Fish not assigned to a PBT population must 
#' be specified as \code{Unassigned}. 
#' 
#' @return a vector
#' 
#' @export

origin_assignment <- function(Clip = indiv$Clip, 
                              PopPBT = indiv$PopPBT) {
  origin = dplyr::if_else(PopPBT == "Unassigned" & Clip == "AI", "wild", 
                   dplyr::if_else(PopPBT == "Unassigned" & Clip == "AD", "hatchery_unknown", "hatchery"))
  origin
} 
