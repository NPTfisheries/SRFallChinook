#functions for fall Chinook run reconstruction

#rate assignment adds the LGR trap rate, the genotyping rate, and the PBT tag rate to the df of individual data
#indiv is a dataframe with individual-level data
#trap_rate_df is a dataframe with the LGR trap rate, genotype rate, and the first date of each new rate period
#pbt_tag_rate is a df with PBT tag rates by pbt population
#sample date is the column in indiv that includes the date when the fish was sampled
#trap_rate_date is the column in trap_rate_df that shows the first date of the new rate period
#trap_rate is the column in trap_rate_df with the LGR trap rate
#geno_rate is the column in trap_rate_df with the genotype sampling rate
#pbt is the column in both indiv and pbt_tag_rate indicating the pbt pop

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
    arrange(StartDate)
  
  indiv %>%
    mutate(TrapRate = trap_rate_df[[trap_rate]][findInterval(.data[[sample_date]], trap_rate_df[[trap_rate_date]])]) %>%
    mutate(GenRate = trap_rate_df[[geno_rate]][findInterval(.data[[sample_date]], trap_rate_df[[trap_rate_date]])]) %>%
    left_join(pbt_tag_rate, by = pbt) 
  
}
  

#length_bin creates length bins to pool fish by fork length
#length bins are in 50mm increments and includes a break point of 570mm to classify jacks in later steps
length_bin <- function(FL = indiv$FL) {
  FL_seq <- seq(20, max(FL, na.rm = T) + 50, by = 50)
  
  pool_FL_interval = cut(FL, breaks = FL_seq, right = FALSE)
  pool_FL = substr(pool_FL_interval, 2, 4)
  #pool_FL = as.numeric(str_extract(pool_FL_interval, "(?<=,)[^)\\]]+"))
  
  pool_FL
}  
  

#expands the count of individual fish to account for the trap rate at LGR
#divides the count of fish by the trap rate

lgr_tr_expansion <- function(count = indiv$count, 
                             trap_rate = indiv$TrapRate) {
  lgr_exp = count/trap_rate
  lgr_exp
} 


#expands the count of individual fish to account for individuals that failed to genotype
#expansion is done proportionally by sex, FL, and trap rate in this example
#indiv is the individual-level dataset
#group_var are the grouping variables; the expansion can be done using different grouping variables
#ng_col is the column containing the information on whether the fish genotyped
#a failed genotype must be coded as "NG"
#count_col is the column containing the counts for individual fish (may have already been expanded using LGR trap rate)

ng_expansion <- function(indiv, 
                         group_var = c("pool_FL", "Sex", "TrapRate", "RelSite"), 
                         ng_col = "RelSite",
                         count_col = lgr_expand) {
  
  pool_var <- setdiff(group_var, ng_col)
  
  #calculate the rate of failed genotyping by length, sex, and trap rate
  ng_prop_tbl <- indiv %>%
    group_by(across(all_of(group_var))) %>%
    summarise(count = sum(.data[[count_col]]), .groups = "drop_last") %>%
    mutate(gen_prop = 1 - (sum(count[.data[[ng_col]] == "NG"])/sum(count))) %>%
    ungroup() %>%
    select(all_of(pool_var), gen_prop) %>%
    distinct()
    
  #expand by the genotyping rate, proportional based on length, sex, and trap rate
  #set NG fish count to 0
  #remove NG fish from the indiv dataset
  ng_expand <- indiv %>%
    full_join(ng_prop_tbl, by = pool_var) %>%
    mutate(ng_expand = .data[[count_col]]/gen_prop) %>%
    mutate(ng_expand = if_else(.data[[ng_col]] == "NG", 0, ng_expand)) %>%
    filter(.data[[ng_col]] != "NG")
  
  ng_prop_tbl
  ng_expand
  
}  


#expand by the PBT tagging rate
#indiv is individual level dataset
#PbtTagRate column containing PBT tag rates
#count_col column containing individual-level counts (may have already been expanded by trap rate, etc.)
#PopPBT column containing assignment to PBT populations. Unassigned must be coded as "Unassigned"
#Clip column specifying ad-intact or ad-clipped. 

pbt_expansion <- function(indiv, 
                          PbtTagRate = "PbtTagRate", 
                          count_col = "ng_expand", 
                          PopPBT = "PopPBT", 
                          Clip = "Clip") {
  
  pbt_expand_tbl <- indiv %>%
    #set PBT tag rate for unassigned fish to 1
    mutate(PbtTagRate = if_else(.data[[PopPBT]] == "Unassigned", 1, .data[[PbtTagRate]])) %>%
    #divide count by tag rate
    mutate(pbt_expand = .data[[count_col]]/PbtTagRate) %>%
    group_by(.data[[Clip]], .data[[PopPBT]]) %>%
    summarise(pbt_expand = sum(pbt_expand), 
              previous_expand = sum(.data[[count_col]]), 
              diff = previous_expand - pbt_expand, .groups = "drop_last") %>%
    group_by(Clip) %>%
    mutate(exp_sum = sum(diff, na.rm = T), 
           diff = if_else(PopPBT == "Unassigned", exp_sum * -1, diff)) %>%
    ungroup() %>%
    mutate(prop_exp = diff/previous_expand)
  
  pbt_expand <- indiv %>%
    left_join(pbt_expand_tbl %>%
                select(Clip, PopPBT, prop_exp), by = c("PopPBT", "Clip")) %>%
    mutate(pbt_expand = .data[[count_col]] - prop_exp) %>%
    select(-prop_exp)
  
  pbt_expand
}

#expand by the genotype sampling rate at LGR
#indiv individual-level dataset
#count_col column with individual-level counts (may have been expanded already)
#gen_rate column with the genotype sampling rate, applied using the rate_assignment function
gen_rate_expansion <- function(indiv, 
                               count_col = "pbt_expand", 
                               gen_rate = "GenRate") {
  gen_rate_expand <- indiv %>%
    mutate(gen_expand = .data[[count_col]]/.data[[gen_rate]])
  
  gen_rate_expand
}

#assign fish as jacks or adults
#adults are >= 570 mm, jacks are < 570 mm
jack_assignment <- function(FL = indiv$pool_FL) {
  jack = if_else(FL >= 570, "adult", "jack")
  jack
}

#assign origin: wild, hatchery assigned, hatchery unknown
#Clip must be coded as "AI" for ad-intact and "AD" for ad-clipped
#PopPBT must be coded as "Unassigned" for unassigned 
origin_assignment <- function(Clip = indiv$Clip, 
                              PopPBT = indiv$PopPBT) {
  origin = if_else(PopPBT == "Unassigned" & Clip == "AI", "wild", 
                   if_else(PopPBT == "Unassigned" & Clip == "AD", "hatchery_unknown", "hatchery"))
  origin
} 


#apply fallback and reascension rates
#fallback rates are assumed to be 0 for wild and hatchery unknown fish
#separate rates for jacks and adults, fish released above and below LGR as juveniles
#fallback is df with fallback_rate defined for jacks and adults released above and below LGR
#jack column specifying jack or adult (must match coding in fallback df)
#origin column specifying origin, hatchery origin must be coded as "hatchery"
#rel_lgr column specifying juvenile release above or below LGR 9must match coding in fallback df)
apply_fallback <- function(indiv, 
                           fallback, 
                           jack = "jack", 
                           origin = "origin", 
                           rel_lgr = "RelLGR",
                           count_col = "gen_expand") {
  
  join_var <- setNames(c(jack, rel_lgr), c(jack, rel_lgr))
  
  fallback_account <- indiv %>%
    left_join(fallback, by = join_var) %>%
    mutate(fallback_rate = as.numeric(fallback_rate)) %>%
    mutate(fallback_rate = if_else(.data[[origin]] != "hatchery", 0, fallback_rate)) %>%
    mutate(fallback_expand = .data[[count_col]]*(1-fallback_rate))
  
  fallback_account
  
}


#age-length key from scale data
#age is df with scale ages
#sex is set to M or F
#sex_col column containing Sex
#scale_age column containing scale ages
#pool_FL column containing binned FL, set using length_bin function

create_age_length_key <- function(age, 
         sex = "M", 
         sex_col = "PhenotypicSex",
         age_col = "Age", 
         pool_FL = "pool_FL") {
  
  age_length_key <- age %>%
    filter(.data[[sex_col]] == sex) %>%
    group_by(.data[[age_col]], .data[[pool_FL]]) %>%
    count() %>% 
    ungroup() %>%
    complete(.data[[age_col]], .data[[pool_FL]], fill = list(n = 0)) %>%
    group_by(.data[[pool_FL]]) %>%
    mutate(prop = n/sum(n)) %>%
    ungroup() %>%
    select(all_of(c(age_col, pool_FL)), prop) %>%
    pivot_wider(names_from = all_of(age_col), values_from = prop) 
  
  age_length_key
}

#assign ages to wild fish
assign_wild_age <- function(wild,
                            age_length_key, 
                            sex = "F",
                            pool_FL = "pool_FL", 
                            count = "count") {
  scale_ages <- colnames(age_length_key)[2:dim(age_length_key)[2]]
  
  wild %>%
    filter(Sex == sex) %>%
    left_join(age_length_key, by = "pool_FL") %>%
    pivot_longer(cols = all_of(scale_ages), names_to = "scale_age", values_to = "prop") %>%
    mutate(count_by_age = prop * .data[[count]]) %>%
    mutate(scale_age = if_else(scale_age == "1", "1.0", scale_age)) %>%
    mutate(age = as.numeric(substr(scale_age, 1, 1)) + as.numeric(substr(scale_age, 3, 3)) + 1) %>%
    mutate(rearing = if_else(substr(scale_age, 1, 1) == 0, "subyearling", "yearling")) %>%
    select(-count, -prop, -scale_age, -agePBT) %>%
    rename(count = count_by_age) %>%
    mutate(brood_year = run_yr - age)
}


#assign ages to hatchery unknown fish
assign_hat_unk_age <- function(hat_unk, 
                               age_length_key, 
                               sex = "F",
                               pool_FL = "pool_FL", 
                               count = "count") {
  
  ages <- colnames(age_length_key)[2:dim(age_length_key)[2]]
  
  hat_unk %>%
    filter(Sex == sex) %>%
    left_join(age_length_key, by = pool_FL) %>%
    pivot_longer(cols = all_of(ages), names_to = "age", values_to = "prop") %>%
    mutate(count_by_age = prop * .data[[count]]) %>%
    mutate(rearing = "subyearling") %>%
    select(-count, -prop, -agePBT)  %>%
    rename(count = count_by_age) %>%
    mutate(age = as.numeric(age)) %>%
    mutate(brood_year = run_yr - age)
}


#trap outage expansion
#beg_outage is the first date of the trap outage, formatted "mm/dd/yyyy"
#end_outage is the last date of the trap outage, formatted "mm/dd/yyyy"
#window count is the total window count during the outage
#final portion of the function accounts for fallback 

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
    filter(DateSampled >= before_outage & DateSampled <= after_outage) %>%
    group_by(.data[[agePBT]], .data[[rel_site]]) %>%
    summarise(count = sum(.data[[count_col]]), .groups = "drop_last") %>%
    ungroup() %>%
    mutate(prop = count/sum(count)) %>%
    mutate(outage_est = prop * window_count) %>%
    select(!count)
  
  #sex proportions from entire run
  sex_ratio <- indiv %>%
    mutate(Sex = if_else(.data[[jack_col]] == "jack", "jack", .data[[sex_col]])) %>%
    group_by(.data[[agePBT]], .data[[rel_site]], .data[[sex_col]]) %>%
    summarise(count = sum(fallback_expand), .groups = "drop_last") %>%
    mutate(sex_prop = count/sum(count)) %>%
    select(!count)
  
  trap_outage_count <- outage_count %>%
    left_join(sex_ratio, by = c(agePBT, rel_site)) %>%
    mutate(outage_count = outage_est * sex_prop) %>%
    select(all_of(c(agePBT, rel_site)), Sex, outage_count)
  
  release_group_locations <- indiv %>%
    select(.data[[rel_site]], RelLGR) %>% distinct()
  
  trap_outage_count <- trap_outage_count %>%
    mutate(origin = if_else(.data[[rel_site]] == "Unassigned", "wild", "hatchery")) %>%
    mutate(jack = if_else(Sex == "jack", "jack", "adult")) %>%
    left_join(release_group_locations, by = rel_site) %>%
    apply_fallback(fallback = fallback, 
                   count_col = "outage_count") %>%
    select(all_of(c(agePBT, rel_site, sex_col)), fallback_expand) %>%
    rename(outage_count = fallback_expand) 
  
}


end_season_expansion <- function(trap_end,
                                  window_count, 
                                  indiv, 
                                  count_col = "pbt_expand",
                                  agePBT = "agePBT", 
                                  rel_site = "RelSite", 
                                  sex_col = "Sex", 
                                  jack_col = "jack") {

  
  before_season_end <- as.Date(trap_end, format = "%m/%d/%Y") - 10
  
  #release group and brood year proportions estimates from period immediately before and after outage
  outage_count <- indiv %>%
    filter(DateSampled >= before_season_end) %>%
    group_by(.data[[agePBT]], .data[[rel_site]]) %>%
    summarise(count = sum(.data[[count_col]]), .groups = "drop_last") %>%  
    ungroup() %>%
    mutate(prop = count/sum(count)) %>%
    mutate(outage_est = prop * window_count) %>%
    select(!count)
  
  #sex proportions from entire run
  sex_ratio <- indiv %>%
    mutate(Sex = if_else(.data[[jack_col]] == "jack", "jack", .data[[sex_col]])) %>%
    group_by(.data[[agePBT]], .data[[rel_site]], .data[[sex_col]]) %>%
    summarise(count = sum(fallback_expand), .groups = "drop_last") %>%
    mutate(sex_prop = count/sum(count)) %>%
    select(!count)
  
  trap_outage_count <- outage_count %>%
    left_join(sex_ratio, by = c(agePBT, rel_site)) %>%
    mutate(outage_count = outage_est * sex_prop) %>%
    select(all_of(c(agePBT, rel_site)), Sex, outage_count)
  
  release_group_locations <- indiv %>%
    select(.data[[rel_site]], RelLGR) %>% distinct()
  
  trap_outage_count <- trap_outage_count %>%
    mutate(origin = if_else(.data[[rel_site]] == "Unassigned", "wild", "hatchery")) %>%
    mutate(jack = if_else(Sex == "jack", "jack", "adult")) %>%
    left_join(release_group_locations, by = rel_site) %>%
    apply_fallback(fallback = fallback, 
                   count_col = "outage_count") %>%
    select(all_of(c(agePBT, rel_site, sex_col)), fallback_expand) %>%
    rename(end_season_count = fallback_expand)         
  
}
