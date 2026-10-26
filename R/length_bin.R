#' @title Creates length bins
#' 
#' @description
#' Creates length bins to pool fish by fork length. Creates 50-mm length bins and includes a break point of 570 mm to classify adults and jacks in later steps. 
#' 
#' @param FL A column containing fork lengths
#' 
#' @return a vector
#' 
#' @export

length_bin <- function(FL = indiv$FL) {
  FL_seq <- seq(20, max(FL, na.rm = T) + 50, by = 50)
  
  pool_FL_interval = cut(FL, breaks = FL_seq, right = FALSE)
  pool_FL = substr(pool_FL_interval, 2, 4)
  
  pool_FL
}  
