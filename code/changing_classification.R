# examine which interactions are classified differently in different islands
df_categorized_check <- distinct(df_categorized)
df_categorized_check <- df_categorized_check %>% select(emln_id, train_layer, test_layer, prop_ones_removed,
                                                        amount_of_removed_1, amount_of_removed_0, prop_0_removed,
                                                        k, lambda, original_links, predicted_values, node_to,
                                                        node_from, removed, input_lambda, itr, predicted_prob_sigm,
                                                        predicted_bin_sigm, original_binary, interaction_id)
identical(unique(df_categorized_check), unique(df_self))

# Identify interactions with inconsistent categories across islands

inconsistent_links <- df_categorized %>%
  group_by(interaction_id) %>%
  summarise(
    n_layers = n_distinct(test_layer),
    n_categories = n_distinct(link_category),
    .groups = "drop"
  ) %>%
  filter(n_categories > 1)

links_with_disagreement <- df_categorized %>%
  semi_join(inconsistent_links, by = "interaction_id")

view(links_with_disagreement %>%
  select(interaction_id, test_layer, link_category) %>%
  distinct() %>%
  arrange(interaction_id, test_layer))

view(distinct(links_with_disagreement)) # in some cases we get different 
# classification for the same interaction in the same island. 
# so we need to average the predicted values of the iterations in df_categorized and only then binarize them
# what is okay: in some cases a locally unique interaction in one island is a possibly missing link in another


# checking if for the same island we can get different classifications
inconsistent_links <- df_categorized_island %>%
  group_by(interaction_id) %>%
  summarise(
    n_layers = n_distinct(test_layer),
    n_categories = n_distinct(link_category),
    .groups = "drop"
  ) %>%
  filter(n_categories > 1)

links_with_disagreement <- df_categorized_island %>%
  semi_join(inconsistent_links, by = "interaction_id")

view(links_with_disagreement %>%
       select(interaction_id, test_layer, link_category) %>%
       distinct() %>%
       arrange(interaction_id, test_layer))

view(distinct(links_with_disagreement)) # in some cases we get different classes

# trying to fix it
df_self_flagged <- df_self %>%
  mutate(
    original_binary = as.integer(original_binary),
    interaction_id  = paste0(node_from, " -> ", node_to)
  ) %>%
  distinct(test_layer, interaction_id, original_binary) %>%
  group_by(interaction_id) %>%
  summarise(
    n_obs_total = sum(original_binary == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  right_join(
    df_self %>%
      mutate(interaction_id = paste0(node_from, " -> ", node_to)) %>%
      distinct(test_layer, interaction_id),
    by = "interaction_id"
  ) %>%
  mutate(
    is_all_zero   = n_obs_total == 0,
    is_unique     = n_obs_total == 1,
    is_shared     = n_obs_total >= 2,
    obs_elsewhere = n_obs_total >= 2  # observed in ≥1 other than focal — approximated
  )

df_removed_flagged <- df_removed %>%
  filter(test_layer == train_layer) %>%
  mutate(
    original_binary    = as.integer(original_binary),
    predicted_bin_sigm = as.integer(predicted_bin_sigm),
    interaction_id     = paste0(node_from, " -> ", node_to)
  ) %>%
  left_join(df_self_flagged, by = c("test_layer", "interaction_id"))
