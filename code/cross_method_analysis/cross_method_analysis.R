# ---- cross method analysis ----
# here we corroborate predictions produced by Barry's application using data from the 
# same plant-frugivore network in Spain that was samples using two different 
# methods.

library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggalluvial)
library(scales)

# remove 20% of the same links in two methods, maybe repeat 50 times 

setwd("~/Documents/github/LP_validation/results/cross_method_results")

# read the two prediction tables. i tried 2 versions: predicting with 20% and 40% link withholding.
mist_nets <- read_csv("hr_mn_40per_prediction_results.csv", show_col_types = FALSE) %>%
  mutate(method = "mist_nets") # m1

observations <- read_csv("hr_obs_40per_prediction_results.csv", show_col_types = FALSE) %>%
  mutate(method = "observations") # m2

## ---- finding missing links using cross-method corroboration ----
# find interactions that are FP in method1 but TP in method2
# the problem is that it's not consistent
fp_to_tp <- mist_nets %>%
  select(higher_level, lower_level, confusion, probability, classification, ground_truth, link) %>%
  rename(confusion_m1 = confusion,
         probability_m1 = probability,
         classification_m1 = classification,
         ground_truth_m1 = ground_truth,
         link_m1 = link) %>%
  inner_join(
    observations %>%
      select(higher_level, lower_level, confusion, probability, classification, ground_truth, link) %>%
      rename(confusion_m2 = confusion,
             probability_m2 = probability,
             classification_m2 = classification,
             ground_truth_m2 = ground_truth,
             link_m2 = link),
    by = c("higher_level", "lower_level")
  ) %>%
  filter(confusion_m1 == "FP", confusion_m2 == "TP")

fp_to_tp # 2 missing links found!
nrow(fp_to_tp)
# write_csv(fp_to_tp, "FP_in_mistnets_TP_in_observations_40per.csv")

# now vise versa
fp_to_tp2 <- observations %>%
  select(higher_level, lower_level, confusion, probability, classification, ground_truth, link) %>%
  rename(confusion_m2 = confusion,
         probability_m2 = probability,
         classification_m2 = classification,
         ground_truth_m2 = ground_truth,
         link_m2 = link) %>%
  inner_join(
    mist_nets %>%
      select(higher_level, lower_level, confusion, probability, classification, ground_truth, link) %>%
      rename(confusion_m1 = confusion,
             probability_m1 = probability,
             classification_m1 = classification,
             ground_truth_m1 = ground_truth,
             link_m1 = link),
    by = c("higher_level", "lower_level")
  ) %>%
  filter(confusion_m2 == "FP", confusion_m1 == "TP")

fp_to_tp2
nrow(fp_to_tp2) # additional missing links found!
# write_csv(fp_to_tp2, "FP_in_observations_TP_in_mistnets_40per.csv")

## ---- full alluvial analysis ----
# withhold the same 20% of links from the two networks
set.seed(42)

# # read matrices
# mistnets_matrix <- read_csv("hr_mn.csv", show_col_types = FALSE)
# observations_matrix <- read_csv("hr_obs.csv", show_col_types = FALSE)
# 
# # check their species are identical
# identical(rownames(mistnets_matrix), rownames(observations_matrix))
# identical(colnames(mistnets_matrix), colnames(observations_matrix))
# dim(mistnets_matrix)
# dim(observations_matrix)
# 
# # keep row names
# rownames(mistnets_matrix) <- mistnets_matrix$AASPECIES
# rownames(observations_matrix) <- observations_matrix$AASPECIES
# 
# mistnets_matrix <- mistnets_matrix %>% select(-AASPECIES)
# observations_matrix <- observations_matrix %>% select(-AASPECIES)
# 
# # find all existing links (non-zero)
# # logical matrix: TRUE only where both have a link
# common_links <- (mistnets_matrix > 0) & (observations_matrix > 0)
# 
# # convert to long format
# common_links_long <- as.data.frame(as.table(common_links)) %>%
#   filter(Freq) %>%
#   rename(
#     lower_level = Var1,
#     higher_level = Var2
#   )
# 
# # randomly select 20% of shared links
# set.seed(42)
# 
# n_remove <- ceiling(0.2 * nrow(common_links_long))
# 
# links_removed <- common_links_long %>%
#   slice_sample(n = n_remove)
# 
# # remove same links from both matrices
# mistnets_removed <- mistnets_matrix
# observations_removed <- observations_matrix
# 
# for (i in seq_len(nrow(links_removed))) {
#   mistnets_removed[
#     links_removed$lower_level[i],
#     links_removed$higher_level[i]
#   ] <- 0
#   
#   observations_removed[
#     links_removed$lower_level[i],
#     links_removed$higher_level[i]
#   ] <- 0
# }

# no link withholding
mist_nets_predictions <- read_csv("hr_mn_no_withholding_prediction_results.csv", show_col_types = FALSE)
observation_predictions <- read_csv("hr_obs_no_withholding_prediction_results.csv", show_col_types = FALSE)

mist_nets_predictions <- mist_nets_predictions %>% mutate(method = "method1")
observation_predictions <- observation_predictions %>% mutate(method = "method2")

df_all <- bind_rows(mist_nets_predictions, observation_predictions) %>%
  mutate(
    interaction_id = paste(higher_level, lower_level, sep = "___")
  )

# compute “observation across methods” flags
obs_counts <- df_all %>%
  group_by(interaction_id) %>%
  summarise(
    n_obs_total = sum(ground_truth == 1),
    .groups = "drop"
  )

# flag links per method
df_flagged <- df_all %>%
  left_join(obs_counts, by = "interaction_id") %>%
  mutate(
    original_binary = ground_truth,
    predicted_bin   = classification,
    
    is_all_zero   = n_obs_total == 0,
    is_unique     = n_obs_total == 1,
    is_shared     = n_obs_total == 2,
    obs_elsewhere = (n_obs_total - (original_binary == 1)) >= 1
  )

# classify links
df_categorized <- df_flagged %>%
  mutate(
    link_category = case_when(
      # observed in only this method
      is_unique   & original_binary == 1 & predicted_bin == 1 ~ "locally_unique_links",
      is_unique   & original_binary == 1 & predicted_bin == 0 ~ "unsupported_links",
      
      # observed in both methods
      is_shared   & original_binary == 1 & predicted_bin == 1 ~ "confirmed_links",
      is_shared   & original_binary == 1 & predicted_bin == 0 ~ "cryptic_links",
      
      # observed nowhere
      is_all_zero & original_binary == 0 & predicted_bin == 0 ~ "likely_forbidden",
      is_all_zero & original_binary == 0 & predicted_bin == 1 ~ "spurious_links",
      
      # observed in the other method only
      obs_elsewhere & original_binary == 0 & predicted_bin == 0 ~ "feasible_links",
      obs_elsewhere & original_binary == 0 & predicted_bin == 1 ~ "possibly_missing_links",
      
      TRUE ~ "unclassified"
    )
  ) 

# summary table
final_table_summary <- df_categorized %>%
  group_by(Category) %>%
  summarise(
    mean = mean(Count), sd = sd(Count),
    min = min(Count), max = max(Count),
    .groups = "drop"
  ) %>%
  arrange(Category)

## ---- plot ----
# make an alluvial plot
# and a heatmap marking where are the missing links
