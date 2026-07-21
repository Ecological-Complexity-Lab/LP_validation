# here we need to start changing. classification should be done per island and per iteration.
# flagging phase is okay because it is done on observations only.
# after that we sum all and calculate percentages for the alluvial plot.
# for the island map, we average the sigmoid probability of the iterations and then classify.
df_self_flagged <- df_self %>%
  left_join(obs_counts, by = "interaction_id") %>%
  mutate(
    is_all_zero   = n_obs_total == 0,
    is_unique     = (n_obs_total == 1),
    is_shared     = (n_obs_total >= 2),
    obs_elsewhere = (n_obs_total - (original_binary == 1)) >= 1
  ) %>%
  select(test_layer, interaction_id, itr, is_all_zero, is_unique, is_shared, obs_elsewhere)


df_removed_flagged <- df_removed %>%
  filter(train_layer == test_layer) %>% 
  mutate(
    original_binary    = as.integer(original_binary),
    predicted_bin_sigm = as.integer(predicted_bin_sigm),
    interaction_id     = paste0(node_from, " -> ", node_to)
  ) %>%
  left_join(df_self_flagged, by = c("test_layer", "interaction_id", "itr")) # now it's okay

# enter the above to main code V
# calculate binary predictions for each island V

# building a data frame for map

# Step 1: Compute average predicted probabilities per interaction per island across iterations
avg_preds <- df_removed_flagged %>%
  group_by(test_layer, interaction_id) %>%
  summarise(
    avg_pred_prob_sigm = mean(predicted_prob_sigm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    avg_bin_sigm = as.integer(avg_pred_prob_sigm > best_discrete_threshold)
  )

# Step 2: Join back to full dataframe and replace predicted_bin_sigm
df_removed_flagged_per_island <- df_removed_flagged %>%
  select(-predicted_bin_sigm) %>%  # drop old version
  left_join(avg_preds, by = c("test_layer", "interaction_id")) %>%
  rename(predicted_bin_sigm = avg_bin_sigm)  # use thresholded avg

df_categorized_per_island <- df_removed_flagged_per_island %>%
  mutate(
    link_category = case_when(
      is_unique   & original_binary == 1 & predicted_bin_sigm == 1 ~ "locally_unique_links",
      is_unique   & original_binary == 1 & predicted_bin_sigm == 0 ~ "unsupported_links",
      is_shared   & original_binary == 1 & predicted_bin_sigm == 1 ~ "confirmed_links",
      is_shared   & original_binary == 1 & predicted_bin_sigm == 0 ~ "cryptic_links",
      is_all_zero & original_binary == 0 & predicted_bin_sigm == 0 ~ "likely_forbidden",
      is_all_zero & original_binary == 0 & predicted_bin_sigm == 1 ~ "spurious_links",
      obs_elsewhere & original_binary == 0 & predicted_bin_sigm == 0 ~ "feasible_links",
      obs_elsewhere & original_binary == 0 & predicted_bin_sigm == 1 ~ "possibly_missing_links",
      TRUE ~ "unclassified"
    )
  )

### ---- identify interactions with inconsistent categories across iterations ----

# Step 1: Detect disagreements across iterations within the same island
inconsistent_links_itr <- df_categorized_per_island %>%
  group_by(test_layer, interaction_id) %>%
  summarise(
    n_iterations = n_distinct(itr),
    n_categories = n_distinct(link_category),
    .groups = "drop"
  ) %>%
  filter(n_categories > 1)

# Step 2: Extract those links with inconsistent classifications
links_with_disagreement_itr <- df_categorized_per_island %>%
  semi_join(inconsistent_links_itr, by = c("test_layer", "interaction_id"))

# Step 3: View inconsistencies
view(
  links_with_disagreement_itr %>%
    select(test_layer, interaction_id, itr, link_category) %>%
    distinct() %>%
    arrange(test_layer, interaction_id, itr)
) # none.

# check disagreements now, in the main table and one island table (notice that you use the average)
# check whether in the alluvial plot we calculate the confusion matrix per island and iteration and sum the confusion matrices

final_table_by_iter_layer <- df_removed_flagged %>%
  group_by(itr, test_layer) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin_sigm == 1),
    FP = sum(original_binary == 0 & predicted_bin_sigm == 1),
    TN = sum(original_binary == 0 & predicted_bin_sigm == 0),
    FN = sum(original_binary == 1 & predicted_bin_sigm == 0),
    
    locally_unique_links = sum(is_unique  & original_binary == 1 & predicted_bin_sigm == 1),
    unsupported_links    = sum(is_unique  & original_binary == 1 & predicted_bin_sigm == 0),
    
    confirmed_links      = sum(is_shared  & original_binary == 1 & predicted_bin_sigm == 1),
    cryptic_links        = sum(is_shared  & original_binary == 1 & predicted_bin_sigm == 0),
    
    likely_forbidden     = sum(is_all_zero & original_binary == 0 & predicted_bin_sigm == 0),
    spurious_links       = sum(is_all_zero & original_binary == 0 & predicted_bin_sigm == 1),
    
    feasible_links       = sum(original_binary == 0 & obs_elsewhere & predicted_bin_sigm == 0),
    possibly_missing_links = sum(original_binary == 0 & obs_elsewhere & predicted_bin_sigm == 1),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = -c(itr, test_layer),
    names_to = "Category", values_to = "Count") %>%
  arrange(test_layer, itr, Category)

final_table_summary <- final_table_by_iter_layer %>%
  group_by(Category) %>%
  summarise(
    total = sum(Count),
    min   = min(Count),
    max   = max(Count),
    sd    = sd(Count),
    .groups = "drop"
  ) %>%
  arrange(Category)


# 1) build flows based on proportions of total sums per island per iteration
cat_totals <- final_table_summary %>%
  transmute(Category, total = coalesce(total, 0)) %>%
  tibble::deframe()

get_total <- function(nm) if (nm %in% names(cat_totals)) unname(cat_totals[[nm]]) else 0

# 2) Build flows: L1 (confusion) → L2 (validation) → L3 (subcategory)
flows <- tibble::tribble(
  ~L1, ~L2,                   ~L3,                     ~value,
  "TP","Validated elsewhere", "confirmed_links",        get_total("confirmed_links"),
  "TP","Not validated",       "locally_unique_links",   get_total("locally_unique_links"),
  "FN","Validated elsewhere", "cryptic_links",          get_total("cryptic_links"),
  "FN","Not validated",       "unsupported_links",      get_total("unsupported_links"),
  "FP","Validated elsewhere", "possibly_missing_links", get_total("possibly_missing_links"),
  "FP","Not validated",       "spurious_links",         get_total("spurious_links"),
  "TN","Validated elsewhere", "feasible_links",         get_total("feasible_links"),
  "TN","Not validated",       "likely_forbidden",       get_total("likely_forbidden")
) %>%
  mutate(
    total_confusion = get_total("TP") + get_total("FP") + get_total("TN") + get_total("FN"),
    prop = ifelse(total_confusion > 0, value / total_confusion, 0),
    alluvium_id = paste(L1, L3, sep = "⟂")
  )



