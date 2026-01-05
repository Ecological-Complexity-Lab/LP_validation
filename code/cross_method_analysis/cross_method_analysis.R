# ---- cross method analysis ----
# here we corroborate predictions produced by Barry's application using data from the 
# same plant-frugivore network in Spain that was samples using two different 
# methods.

library(dplyr)
library(readr)

# read the two prediction tables
mist_nets <- read_csv("hr_mn_prediction_results.csv", show_col_types = FALSE) %>%
  mutate(method = "mist_nets") # m1

observations <- read_csv("hr_obs_prediction_results.csv", show_col_types = FALSE) %>%
  mutate(method = "observations") # m2

# find interactions that are FP in method1 but TP in method2
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

fp_to_tp
# write_csv(fp_to_tp, "FP_in_method1_TP_in_method2.csv")
