#' @title Assign trap rates, genotype rates, and PBT tag rates to individual data
#' 
#' @description
#' Adds the Lower Granite trap rates, genotype sample rates, and PBT tagging rates to the dataframe with the individual-level data from the Lower Granite Dam adult trap. 
#' 
#' @param indiv A dataframe with the individual-level data from the Lower Granite Dam adult trap. This dataframe must contain a column with the date the fish was sampled 
#' and a column with the population it was assigned to using PBT. 
#' 
#' @param trap_rate_df A dataframe containing data on the operations of the Lower Granite Dam adult trap. This dataframe must contain a column with the trapping rate, 
#' a column with the first date of each trap rate period, and the genotyping rate for each trapping period. 
#' 
#' @param pbt_tag_rate A dataframe containing each PBT population and the associated tag rate.
#' 
#' @param sample_date The column in \code{indiv} with the date each fish was sampled.
#' 
#' @param trap_rate_date The column in \code{trap_rate_df} with the first date of each trap rate period.
#' 
#' @param trap_rate The column in \code{trap_rate_df} with the trap rate for each period.
#' 
#' @param geno_rate The column in \code{trap_rate_df} with the genotyping rate for each period. 
#' 
#' @param pbt The name of the column, shared by \code{indiv} and \code{pbt_tag_rate}, containing the PBT population assignments.
#' 
#' @return A dataframe
#' 
#' @export

rate_assignment <- function(indiv, 
                            trap_rate_df, 
                            pbt_tag_rate, 
                            sample_date = "DateSampled", 
                            trap_rate_date = "StartDate", 
                            trap_rate = "TrapRate", 
                            geno_rate = "geno_rate", 
                            pbt = "PopPBT")
{
  #arrange trap rate by date
  trap_rate_df <- trap_rate_df %>%
    dplyr::arrange(StartDate)
  
  indiv %>%
    dplyr::mutate(TrapRate = trap_rate_df[[trap_rate]][findInterval(.data[[sample_date]], trap_rate_df[[trap_rate_date]])]) %>%
    dplyr::mutate(GenRate = trap_rate_df[[geno_rate]][findInterval(.data[[sample_date]], trap_rate_df[[trap_rate_date]])]) %>%
    dplyr::left_join(pbt_tag_rate, by = pbt) 
  
}
