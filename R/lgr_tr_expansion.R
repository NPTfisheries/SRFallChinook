#' @title Apply trap rate expansion
#' 
#' @description
#' Applies the expansion for the trap rate at Lower Granite Dam by dividing the count of fish by the trap rate. Typically, this expansion is applied first, so each row represents one fish, 
#' and each individual fish is expanded by the appropriate trap rate. 
#' 
#' @param count A column containing a count for each row. If this expansion is applied first, \code{count} will be \code{1} for each fish. 
#' 
#' @param trap_rate A column containing the corresponding trap rate at the Lower Granite Dam adult trap. 
#' 
#' @return a vector
#' 
#' @export


lgr_tr_expansion <- function(count = indiv$count, 
                             trap_rate = indiv$TrapRate) {
  lgr_exp = count/trap_rate
  lgr_exp
} 
