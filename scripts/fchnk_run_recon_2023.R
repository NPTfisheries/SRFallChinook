
indiv <- indiv %>%
  mutate(PopPBT = if_else(PopPBT == "OtsPBt19-PRHFa", "OtsPBT19-PRHFa", PopPBT))

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

#2023 does not have raw age data

#use scale ages to create age-length keys for wild fish
# age <- age %>%
#   mutate(pool_FL = length_bin(FL*10), 
#          pool_FL = as.numeric(pool_FL), 
#          #to make length bins match spreadsheet
#          pool_FL = if_else(pool_FL == FL*10, pool_FL, pool_FL + 50))
# 
# wild_age_len_key <- c("M", "F") %>%
#   set_names() %>%
#   map(~ create_age_length_key(age,
#                               sex = .x,
#                               sex_col = "GenSex",
#                               age_col = "Age",
#                               pool_FL = "pool_FL"))

#read in age-length key directly and format to match function output
f_wild_alk <- read.csv("data_inputs/2023/wild_female_age_length_23.csv", check.names = FALSE)
m_wild_alk <- read.csv("data_inputs/2023/wild_male_age_length_23.csv", check.names = FALSE)

names(f_wild_alk)[1] <- "pool_FL"
names(m_wild_alk)[1] <- "pool_FL"

wild_age_len_key <- list("M" = m_wild_alk, "F" = f_wild_alk)



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
  
#No outages in 2023, trapping continued until 11/15/23

#account for trap outage expansions
#mid-season outage
# beg_outage <- "09/16/2025"
# end_outage <- "09/20/2025"
# window_count <- 3230
# 
# trap_outage_count <- trap_outage_expansion(beg_outage = beg_outage, 
#                                            end_out = end_outage, 
#                                            window_count = window_count, 
#                                            indiv = indiv, 
#                                            count_col = "pbt_expand", 
#                                            agePBT = "agePBT", 
#                                            rel_site = "RelSite", 
#                                            sex_col = "Sex", 
#                                            jack_col = "jack")
# 
# #trapping ends before the end of the run
# end_season_count <- end_season_expansion(trap_end = "11/04/2025", 
#                                              window_count = 1064, 
#                                              indiv = indiv, 
#                                              count_col = "pbt_expand", 
#                                              agePBT = "agePBT", 
#                                              rel_site = "RelSite", 
#                                              sex_col = "Sex", 
#                                              jack_col = "jack")


#hatchery assigned fish can be added to hat_by_rel_site df
# hat_by_rel_site_with_outage <- hat_by_rel_site %>%
#   rename(agePBT = age) %>%
#   left_join(trap_outage_count, by = c("agePBT", "RelSite", "Sex")) %>%
#   left_join(end_season_count, by = c("agePBT", "RelSite", "Sex")) %>%
#   replace_na(list(outage_count = 0, end_season_count = 0)) %>%
#   mutate(total = count + outage_count + end_season_count) %>%
#   rename(age = agePBT) %>%
#   select(-outage_count, -end_season_count)

#wild fish need to have ages assigned based on total run proportions
#spreadsheet says outage age proportions are applied to wild fish based on time period
#this part is a bit messy and could use some clean up
# wild_with_outage <- wild_by_age %>%
#   mutate(prop = count/sum(count)) %>%
#   left_join(trap_outage_count %>% filter(RelSite == "Unassigned"), by = c("Sex", "RelSite")) %>%
#   left_join(end_season_count %>% filter(RelSite == "Unassigned"), by = c("Sex", "RelSite")) %>%
#   mutate(outage_expand = outage_count + end_season_count) %>%
#   mutate(outage_n_by_age = count*prop) %>%
#   mutate(total = count + outage_n_by_age) %>%
#   select(RelSite, Sex, age, rearing, brood_year, origin, Hatchery, count, total)

hat_by_rel_site <- hat_by_rel_site %>%
  rename(total = count) %>%
  #add rearing info for 2 release groups not in the release_group df
  mutate(rearing = if_else(RelSite %in% c("LFH UNK", "LFHSO2"), "subyearling", rearing)) 

wild_by_age <- wild_by_age %>%
  rename(total = count)

hat_unk_by_age <- hat_unk %>%
  mutate(Sex = if_else(jack == "jack", "jack", Sex)) %>%
  group_by(RelSite, Sex, age, rearing, brood_year, origin, Hatchery) %>%
  summarise(count = sum(count, na.rm = T), .groups = "drop") %>%
  rename(total = count)
  

final_estimates <- rbind(hat_by_rel_site, 
                         wild_by_age, 
                         hat_unk_by_age)

final_estimates <- final_estimates %>%
  mutate(RelSite = if_else(RelSite == "Unassigned AD", "Unassigned", RelSite))

#saveRDS(final_estimates, file = "output/final_estimates.RDS")

final_totals <- read.csv("data_inputs/2023/final_totals_check_23.csv")
final_totals <- final_totals %>%
  rename('F' = Female, 
         M = Males, 
         jack = jacks) %>%
  pivot_longer(cols = c("F", "M", "jack"), names_to = "Sex", values_to = "spreadsheet_count") %>%
  rename(age = agePBT, 
         RelSite = grp, 
         brood_year = Brd.Yr) %>%
  mutate(RelSite = if_else(RelSite == "Wild", "Unassigned", 
                           if_else(RelSite == "Unkn Hat AD", "Unassigned", RelSite))) %>%
  select(-totals)


check_final <- full_join(final_estimates, final_totals, by = c("Sex", "RelSite", "brood_year", "age", "rearing" , "origin"))

check_final <- check_final %>%
  mutate(total = if_else(is.na(total), 0, total)) %>%
  mutate(spreadsheet_count = if_else(is.na(spreadsheet_count), 0, spreadsheet_count)) %>%
  mutate(diff = total - spreadsheet_count) %>%
  mutate(pct_diff = (diff/spreadsheet_count)*100)

summary(check_final$diff)

  

raw_diff_hist <- ggplot(check_final, aes(x = diff)) + 
  geom_histogram(fill = "grey", color = "black", binwidth = 5)  + 
  plot_theme + 
  xlab("Difference between estimates") + 
  ylab("Count")  
  #scale_x_continuous(expand = c(0, 0), limits = c(-50, 150))
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 125))


final_total_plot <- ggplot(check_final, aes(x = spreadsheet_count, y = total, group = origin)) + 
  geom_abline(slope = 1, linewidth = 1) +
  geom_point(size = 4, shape = 21, aes(fill = origin), color = "black") + 
  plot_theme + 
  xlab("Original method") + 
  ylab("New method") 
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 2000)) + 
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 2000))

check_final %>%
  filter(origin == "hatchery") %>%
  ggplot(aes(x = spreadsheet_count, y = total)) + 
  geom_point()

#total escapement estimates are pretty close, ~1000 off
#age and origin specific estimates are not good
summary_table <- check_final %>%
  group_by(origin) %>%
  summarise('New method' = round(sum(total, na.rm = T)), 'Original method' = sum(spreadsheet_count, na.rm = T))



hist_by_origin <- ggplot(check_final, aes(x = diff, group = origin)) + 
  geom_histogram(aes(fill = origin), color = "black", binwidth = 5)  + 
  plot_theme + 
  xlab("Difference between estimates") + 
  ylab("Count") + 
  #scale_x_continuous(expand = c(0, 0), limits = c(-50, 150)) + 
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + 
  facet_wrap(~origin) + 
  theme(strip.background = element_rect(fill = "grey70", color = "black"), 
        strip.text = element_text(size = 12), 
        panel.spacing = unit(10, "mm")) + 
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dashed")
