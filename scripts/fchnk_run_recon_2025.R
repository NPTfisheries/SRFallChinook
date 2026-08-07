
#assign LGR trap rate, genotype sample rates, and pbt tagging rates to the individual-level data
indiv <- rate_assignment(indiv, 
                trap_rate_df = trap_rate,
                pbt_tag_rate = pbt_tag_rate, 
                sample_date = "DateSampled",
                trap_rate_date = "StartDate", 
                trap_rate = "TrapRate", 
                geno_rate = "geno_rate", 
                pbt = "PopPBT")

#pool length bins
indiv <- indiv %>%
  mutate(pool_FL = length_bin(FL)) 

#this is to get my length bins to match the length bins in the spreadsheet
#except for 990 and 1000 mm length bins
indiv <- indiv %>%
  mutate(pool_FL = as.numeric(pool_FL), 
         pool_FL = if_else(pool_FL == FL, pool_FL, pool_FL + 50))


#expand based on rate of failed genotyping 
indiv <- ng_expansion(indiv, 
             group_var = c("pool_FL", "Sex", "TrapRate", "RelSite"), 
             ng_col = "RelSite", 
             count_col = "count")

#expand based on Lower Granite trap rate
indiv <- indiv %>%
  mutate(lgr_expand = lgr_tr_expansion(ng_expand, TrapRate))



#expand based on PBT tag rate
indiv <- pbt_expansion(indiv, 
                       PbtTagRate = "PbtTagRate", 
                       count_col = "lgr_expand", 
                       PopPBT = "PopPBT", 
                       Clip = "Clip")

#expand based on rate of genotype sampling
indiv <- gen_rate_expansion(indiv, 
                            count_col = "pbt_expand", 
                            gen_rate = "GenRate")

#assign jacks vs adults
indiv <- indiv %>%
  mutate(jack = jack_assignment(pool_FL))

#assign origin (hatchery, wild, hatchery unknown)
indiv <- indiv %>%
  mutate(origin = origin_assignment(Clip, PopPBT))

#account for fallback and reascension
indiv <- apply_fallback(indiv, fallback, jack = "jack", 
               origin = "origin", 
               count_col = "gen_expand", 
               rel_lgr = "RelLGR")

#use scale ages to create age-length keys for wild fish
age <- age %>%
  mutate(pool_FL = length_bin(FL*10), 
         pool_FL = as.numeric(pool_FL), 
         #to make length bins match spreadsheet
         pool_FL = if_else(pool_FL == FL*10, pool_FL, pool_FL + 50))

wild_age_len_key <- c("M", "F") %>%
  set_names() %>%
  map(~ create_age_length_key(age,
                              sex = .x,
                              sex_col = "GenSex",
                              age_col = "Age",
                              pool_FL = "pool_FL"))

#use PBT ages to create age-length key for hatchery fish
#to apply to hatchery unknown
hat_unk_age_len_key <- c("M", "F") %>%
  set_names() %>%
  map(~ create_age_length_key(subset(indiv, origin == "hatchery"), 
                              sex = .x,
                              sex_col = "Sex",
                              age_col = "agePBT",
                              pool_FL = "pool_FL"))


#assign ages for wild and hatchery unknown fish
grouped <- indiv %>% group_by(agePBT, 
                   RelSite, 
                   Sex, 
                   Hatchery, 
                   pool_FL, 
                   jack, 
                   origin) %>%
  summarise(count = sum(fallback_expand, na.rm = T), .groups = "drop_last") %>%
  ungroup()

#wild fish
wild <- grouped %>%
  filter(origin == "wild")

wild <- c("M", "F") %>%
  set_names() %>%
  map_dfr(~ assign_wild_age(
    wild, 
    sex = .x,
    age_length_key = wild_age_len_key[[.x]], 
    pool_FL = "pool_FL", 
    count = "count"
))

#hatchery unknown fish
hat_unk <- grouped %>%
  filter(origin == "hatchery_unknown")

hat_unk <- c("M", "F") %>%
  set_names() %>%
  map_dfr(~ assign_hat_unk_age(
    hat_unk, 
    sex = .x,
    age_length_key = hat_unk_age_len_key[[.x]], 
    pool_FL = "pool_FL", 
    count = "count"
  ))

#hatchery assigned 
hat_assigned <- grouped %>%
  filter(origin == "hatchery") %>%
  rename(age = agePBT) %>%
  mutate(brood_year = run_yr - age) %>%
  left_join(release_group %>%
              select(brood_year, RelSite, rearing), by = c("brood_year", "RelSite")) %>%
  #assume all strays are subyearling releases
  mutate(rearing = if_else(Hatchery %in% SR_hatchery, rearing, "subyearling")) %>%
  #add rearing for BY2023 LFYO (missing from release_group data) 
  mutate(rearing = if_else(RelSite == "LFYO" & brood_year == 2023, "yearling", rearing))

hat_by_rel_site <- hat_assigned %>%
  mutate(Sex = if_else(jack == "jack", "jack", Sex)) %>%
  group_by(RelSite, Sex, age, rearing, brood_year, origin, Hatchery) %>%
  summarise(count = sum(count, na.rm = T), .groups = "drop")

# est_by_FL <- rbind(wild, hat_unk, hat_assigned)
# 
# est_by_rel_site <- est_by_FL %>%
#   mutate(Sex = if_else(jack == "jack", "jack", Sex)) %>%
#   group_by(RelSite, Sex, age, rearing, brood_year, origin) %>%
#   summarise(count = sum(count, na.rm = T), .groups = "drop")

wild_by_age <- wild %>%
  mutate(Sex = if_else(jack == "jack", "jack", Sex)) %>%
  group_by(RelSite, Sex, age, rearing, brood_year, origin, Hatchery) %>%
  summarise(count = sum(count, na.rm = T), .groups = "drop")
  

#account for trap outage expansions
#mid-season outage
beg_outage <- "09/16/2025"
end_outage <- "09/20/2025"
window_count <- 3230

trap_outage_count <- trap_outage_expansion(beg_outage = beg_outage, 
                                           end_out = end_outage, 
                                           window_count = window_count, 
                                           indiv = indiv, 
                                           count_col = "pbt_expand", 
                                           agePBT = "agePBT", 
                                           rel_site = "RelSite", 
                                           sex_col = "Sex", 
                                           jack_col = "jack")

#trapping ends before the end of the run
end_season_count <- end_season_expansion(trap_end = "11/04/2025", 
                                             window_count = 1064, 
                                             indiv = indiv, 
                                             count_col = "pbt_expand", 
                                             agePBT = "agePBT", 
                                             rel_site = "RelSite", 
                                             sex_col = "Sex", 
                                             jack_col = "jack")


#hatchery assigned fish can be added to hat_by_rel_site df
hat_by_rel_site_with_outage <- hat_by_rel_site %>%
  rename(agePBT = age) %>%
  left_join(trap_outage_count, by = c("agePBT", "RelSite", "Sex")) %>%
  left_join(end_season_count, by = c("agePBT", "RelSite", "Sex")) %>%
  replace_na(list(outage_count = 0, end_season_count = 0)) %>%
  mutate(total = count + outage_count + end_season_count) %>%
  rename(age = agePBT) %>%
  select(-outage_count, -end_season_count)

#wild fish need to have ages assigned based on total run proportions
#spreadsheet says outage age proportions are applied to wild fish based on time period
#this part is a bit messy and could use some clean up
wild_with_outage <- wild_by_age %>%
  mutate(prop = count/sum(count)) %>%
  left_join(trap_outage_count %>% filter(RelSite == "Unassigned"), by = c("Sex", "RelSite")) %>%
  left_join(end_season_count %>% filter(RelSite == "Unassigned"), by = c("Sex", "RelSite")) %>%
  mutate(outage_expand = outage_count + end_season_count) %>%
  mutate(outage_n_by_age = count*prop) %>%
  mutate(total = count + outage_n_by_age) %>%
  select(RelSite, Sex, age, rearing, brood_year, origin, Hatchery, count, total)

hat_unk_by_age <- hat_unk %>%
  mutate(Sex = if_else(jack == "jack", "jack", Sex)) %>%
  group_by(RelSite, Sex, age, rearing, brood_year, origin, Hatchery) %>%
  summarise(count = sum(count, na.rm = T), .groups = "drop") %>%
  mutate(total = count)
  

final_estimates <- rbind(hat_by_rel_site_with_outage, 
                         wild_with_outage, 
                         hat_unk_by_age)








#check spreadsheet estimates

#match 'Hat by PBT' sheet
hat_by_pbt_chk <- hat_assigned %>%
  mutate(Sex = if_else(jack == "jack", "jack", Sex)) %>%
  select(!jack) %>%
  group_by(age, RelSite, Sex, brood_year, rearing) %>%
  summarise(total = sum(count))

hat_by_pbt_spreadsheet <- read.csv("data_inputs/hat_by_pbt_check.csv")
hat_by_pbt_spreadsheet <- hat_by_pbt_spreadsheet %>%
  pivot_longer(cols = c("F", "M", "jack"), names_to = "Sex", values_to = "spreadsheet_count")

check_output_hatbypbt <- full_join(hat_by_pbt_chk, hat_by_pbt_spreadsheet, 
                                   by = c("Sex", "brood_year", "RelSite"))  

#output matches for hatchery assigned fish by age, release group, and brood year
ggplot(check_output_hatbypbt, aes(x = total, y = spreadsheet_count)) + 
  geom_point()

#match final totals from spreadsheet
final_totals <- read.csv("data_inputs/final_totals_check.csv")
final_totals <- final_totals %>%
  pivot_longer(cols = c("F", "M", "jack"), names_to = "Sex", values_to = "spreadsheet_count")


check_final <- full_join(final_estimates, final_totals, by = c("Sex", "RelSite", "brood_year", "age", "rearing" , "origin"))

ggplot(check_final, aes(x = spreadsheet_count, y = total)) + 
  geom_point() + 
  geom_abline(slope = 1)

check_final <- check_final %>%
  mutate(diff = total - spreadsheet_count) %>%
  mutate(pct_diff = (diff/spreadsheet_count)*100)

summary(check_final$diff)

ggplot(check_final, aes(x = pct_diff)) + 
  geom_histogram(fill = "grey", color = "black")

ggplot(check_final, aes(x = pct_diff)) + 
  geom_histogram(fill = "grey", color = "black") + 
  facet_wrap(~origin)

#My method seems to have slightly underestimated wild fish abundance
#need to look at where things went wrong
#likely somewhere in age assignment and particularly age assignment for outage expansion fish

check_final %>%
  group_by(RelSite, Sex, rearing, origin) %>%
  summarise(r = sum(total, na.rm = T), excel = sum(spreadsheet_count, na.rm = T)) %>%
  ggplot(aes(x = r, y = excel)) + 
  geom_point() + 
  geom_abline(slope = 1)
