#' @title Classifies fish as adults or jacks
#' 
#' @description
#' Fish with fork lengths of 570mm or less are classified as jacks; fish with fork lengths >570mm are classified as adults.  
#' 
#' @param FL A column containing fork lengths or length bins
#' 
#' @return a vector
#' 
#' @export


jack_assignment <- function(FL = indiv$pool_FL) {
  jack = dplyr::if_else(FL >= 570, "adult", "jack")
  jack
}
