# ---- Serra-Martin link classification and alluvial plots ----
# Loads cross-site SVD prediction results, classifies predicted links via
# spatial corroboration across 6 habitat patches, then validates each method's
# classifications using evidence from the other sampling method.
# Produces two 4-axis alluvial plots (one per method) and a tile map for
# the highest-connectance site, following canary_link_classification.R.

source("code/common.R")
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggalluvial)
library(scales)

## ---- 1. Alluvial aesthetics (matching cross_location_cross_method.R) ----
col_confusion <- c("TP" = "lightsteelblue",  "FP" = "lightsteelblue2",
                   "TN" = "rosybrown",        "FN" = "rosybrown2")
col_validation <- c("Observed elsewhere"     = "sandybrown",
                    "Not observed elsewhere" = "thistle3")
col_subcats <- c(
  "recurrent"        = "coral3",
  "locally_unique"   = "thistle3",
  "model_elusive"    = "coral",
  "unsupported"      = "thistle1",
  "possibly_missing" = "coral2",
  "unconfirmed"      = "thistle",
  "locally_absent"   = "coral1",
  "likely_forbidden" = "thistle2"
)
stratum_fill <- c(col_confusion, col_validation, col_subcats)

order_confusion  <- c("TP", "FP", "TN", "FN")
order_validation <- c("Not observed elsewhere", "Observed elsewhere")
order_subtypes   <- c(
  "locally_unique",   "unconfirmed",
  "likely_forbidden", "unsupported",
  "recurrent",        "possibly_missing",
  "locally_absent",   "model_elusive"
)

metric          <- "prop"
flow_alpha      <- 0.7
flow_colour     <- NA
bg_col          <- "white"
name_size       <- 4.2
label_nudge_x   <- 0.03
label_color_map <- c("Confusion"  = "mistyrose4",
                     "Validation" = "mistyrose4",
                     "Subtype"    = "mistyrose4")

## ---- 2. Functions (from cross_location_cross_method.R) ----
classify_by_location <- function(df_all, method_name) {
  df_method <- df_all %>% filter(method == method_name)

  obs_counts <- df_method %>%
    group_by(interaction_id) %>%
    summarise(n_obs_total = sum(ground_truth == 1), .groups = "drop")

  df_flagged <- df_method %>%
    left_join(obs_counts, by = "interaction_id") %>%
    mutate(
      original_binary = ground_truth,
      predicted_bin   = classification,

      is_all_zero   = n_obs_total == 0,
      is_unique     = n_obs_total == 1,
      is_shared     = n_obs_total >= 2,
      obs_elsewhere = (n_obs_total - as.integer(original_binary == 1)) >= 1,
      validation    = ifelse(obs_elsewhere, "Observed elsewhere", "Not observed elsewhere")
    )

  df_categorized <- df_flagged %>%
    mutate(
      link_category = case_when(
        is_unique     & original_binary == 1 & predicted_bin == 1 ~ "locally_unique",
        is_unique     & original_binary == 1 & predicted_bin == 0 ~ "unsupported",
        is_shared     & original_binary == 1 & predicted_bin == 1 ~ "recurrent",
        is_shared     & original_binary == 1 & predicted_bin == 0 ~ "model_elusive",
        is_all_zero   & original_binary == 0 & predicted_bin == 0 ~ "likely_forbidden",
        is_all_zero   & original_binary == 0 & predicted_bin == 1 ~ "unconfirmed",
        obs_elsewhere & original_binary == 0 & predicted_bin == 0 ~ "locally_absent",
        obs_elsewhere & original_binary == 0 & predicted_bin == 1 ~ "possibly_missing",
        TRUE ~ "unclassified"
      ),
      confusion = case_when(
        original_binary == 1 & predicted_bin == 1 ~ "TP",
        original_binary == 1 & predicted_bin == 0 ~ "FN",
        original_binary == 0 & predicted_bin == 1 ~ "FP",
        original_binary == 0 & predicted_bin == 0 ~ "TN",
        TRUE ~ "UNK"
      )
    )

  summary_tbl <- df_categorized %>%
    filter(link_category != "unclassified", confusion != "UNK") %>%
    count(location, confusion, validation, link_category, name = "value")

  list(categorized = df_categorized, summary = summary_tbl)
}

make_alluvial_validated <- function(df_categorized, add_obs_ids, method_label,
                                    validation_label = "Additional method") {
  order_L4 <- as.vector(rbind(
    paste0(order_subtypes, " — confirmed"),
    paste0(order_subtypes, " — unconfirmed")
  ))
  col_L4 <- c(
    setNames(unname(col_subcats[order_subtypes]),
             paste0(order_subtypes, " — confirmed")),
    setNames(rep("grey88", length(order_subtypes)),
             paste0(order_subtypes, " — unconfirmed"))
  )
  stratum_fill_4    <- c(stratum_fill, col_L4)
  axis4_label       <- validation_label
  order_axis4       <- c("Confusion", "Validation", "Subtype", axis4_label)
  label_color_map_4 <- c(label_color_map, setNames("mistyrose4", axis4_label))

  flows_tbl <- df_categorized %>%
    filter(link_category != "unclassified", confusion != "UNK") %>%
    mutate(
      add_obs = interaction_id %in% add_obs_ids,
      L4      = paste0(link_category,
                       ifelse(add_obs, " — confirmed", " — unconfirmed"))
    ) %>%
    count(confusion, validation, link_category, L4, name = "value") %>%
    mutate(
      L1          = confusion,
      L2          = validation,
      L3          = link_category,
      total       = sum(value),
      prop        = ifelse(total > 0, value / total, 0),
      alluvium_id = paste(L1, L3, L4, sep = "⟂")
    )

  flows_long <- flows_tbl %>%
    select(L1, L2, L3, L4, value, prop, alluvium_id) %>%
    tidyr::pivot_longer(c(L1, L2, L3, L4), names_to = "axis", values_to = "stratum") %>%
    dplyr::mutate(
      axis = dplyr::recode(axis,
        L1 = "Confusion", L2 = "Validation", L3 = "Subtype", L4 = axis4_label
      ),
      stratum = dplyr::case_when(
        axis == "Confusion"  ~ factor(stratum, levels = order_confusion),
        axis == "Validation" ~ factor(stratum, levels = order_validation),
        axis == "Subtype"    ~ factor(stratum, levels = order_subtypes),
        axis == axis4_label  ~ factor(stratum, levels = order_L4),
        TRUE ~ factor(stratum)
      ),
      axis = factor(axis, levels = order_axis4)
    )

  yr <- diff(range(flows_long[[metric]], na.rm = TRUE))

  stratum_totals <- flows_long %>%
    dplyr::group_by(axis, stratum) %>%
    dplyr::summarise(val = sum(.data[[metric]]), .groups = "drop") %>%
    dplyr::mutate(stratum_chr = as.character(stratum))

  tmp_build <- ggplot_build(
    ggplot(flows_long,
           aes(x = axis, stratum = stratum, alluvium = alluvium_id,
               y = .data[[metric]])) +
      geom_stratum(width = 0.03, color = NA)
  )

  geo       <- as.data.frame(tmp_build$data[[1]])
  ax_levels <- levels(flows_long$axis)

  label_geom <- geo %>%
    dplyr::transmute(
      x_mid       = (xmin + xmax) / 2,
      y_mid       = (ymin + ymax) / 2,
      stratum_chr = as.character(stratum),
      axis        = ax_levels[round(x)]
    )

  label_df <- dplyr::left_join(stratum_totals, label_geom,
                               by = c("axis", "stratum_chr")) %>%
    dplyr::mutate(
      name_txt    = stratum_chr %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_to_sentence(),
      pct_txt     = if (metric == "prop") scales::percent(val, accuracy = 1)
                    else                  scales::label_number_si()(val),
      label_final = paste0(name_txt, " (", pct_txt, ")")
    ) %>%
    dplyr::arrange(axis, y_mid) %>%
    dplyr::group_by(axis) %>%
    dplyr::mutate(
      label_y = y_mid + (row_number() - mean(row_number())) * (0.033 * yr)
    ) %>%
    dplyr::ungroup()

  ggplot(
    flows_long,
    aes(x = axis, stratum = stratum, alluvium = alluvium_id,
        y = .data[[metric]], fill = stratum)
  ) +
    geom_stratum(width = 0.02, color = NA, fill = bg_col, alpha = 1) +
    geom_stratum(width = 0.03, color = "white") +
    geom_alluvium(color = flow_colour, alpha = flow_alpha,
                  width = 0.15, knot.pos = 0.2) +
    scale_fill_manual(values = stratum_fill_4, guide = "none") +
    scale_y_continuous(
      labels = if (metric == "prop") percent_format(accuracy = 1)
               else                  label_number_si()
    ) +
    geom_text(
      data        = label_df,
      inherit.aes = FALSE,
      aes(x = x_mid + label_nudge_x, y = label_y,
          label = label_final, color = axis),
      fontface = "bold", size = name_size, hjust = 0
    ) +
    scale_color_manual(values = label_color_map_4, guide = "none") +
    coord_cartesian(clip = "off") +
    labs(
      title    = paste("Link classification + additional method validation —", method_label),
      subtitle = paste("Confusion → Spatial corroboration → Subtype →", validation_label),
      caption  = "Proportions reflect per-location assessments — each interaction is assessed once per site.",
      x        = NULL,
      y        = if (metric == "prop") "Proportion of interactions" else "Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(size = 11, face = "bold"),
      plot.title         = element_text(face = "bold"),
      plot.subtitle      = element_text(color = "grey30"),
      plot.caption       = element_text(color = "grey50", size = 8, hjust = 0),
      plot.margin        = margin(20, 220, 20, 20)
    )
}

## ---- 3. Load prediction CSVs ----
df_all <- read_csv(
  "results/predictions/serra_martin_loo_prediction_results.csv",
  show_col_types = FALSE
) %>%
  rename(location = focal_site, classification = prediction) %>%
  mutate(interaction_id = paste(higher_level, lower_level, sep = "___"))

## ---- 4. Classify by location (spatial corroboration) ----
results_obs <- classify_by_location(df_all, "obs")
results_rpi <- classify_by_location(df_all, "rpi")

df_obs_categorized <- results_obs$categorized
summary_obs        <- results_obs$summary

df_rpi_categorized <- results_rpi$categorized
summary_rpi        <- results_rpi$summary

## ---- 5. Cross-method validation ----
# Additional evidence: interaction IDs with ground_truth == 1 in the other method
add_obs_for_obs <- df_all %>%
  filter(method == "rpi", ground_truth == 1) %>%
  pull(interaction_id) %>% unique()

add_obs_for_rpi <- df_all %>%
  filter(method == "obs", ground_truth == 1) %>%
  pull(interaction_id) %>% unique()

cat(sprintf("Cross-method (rpi → obs validation): %d interaction IDs\n",
            length(add_obs_for_obs)))
cat(sprintf("Cross-method (obs → rpi validation): %d interaction IDs\n",
            length(add_obs_for_rpi)))

## ---- 6. Alluvial plots ----
gg_obs_validated <- make_alluvial_validated(
  df_obs_categorized,
  add_obs_ids      = add_obs_for_obs,
  method_label     = "Direct observation (obs)",
  validation_label = "Raspberry Pi camera (rpi)"
)
gg_obs_validated

gg_rpi_validated <- make_alluvial_validated(
  df_rpi_categorized,
  add_obs_ids      = add_obs_for_rpi,
  method_label     = "Raspberry Pi camera (rpi)",
  validation_label = "Direct observation (obs)"
)
gg_rpi_validated

## ---- 7. Highest-connectance site ----
connectance_table <- df_all %>%
  group_by(location, method) %>%
  summarise(
    n_observed    = sum(ground_truth == 1),
    n_plants      = n_distinct(lower_level),
    n_pollinators = n_distinct(higher_level),
    possible_links = n_plants * n_pollinators,
    connectance   = n_observed / possible_links,
    .groups = "drop"
  ) %>%
  arrange(desc(connectance))

print(connectance_table)

top_row    <- connectance_table %>% slice_max(connectance, n = 1)
top_site   <- top_row$location
top_method <- top_row$method

cat(sprintf("\nHighest-connectance site: %s (%s) — connectance = %.3f\n",
            top_site, top_method, top_row$connectance))

## ---- 8. Map for highest-connectance site ----
df_site_cat <- (if (top_method == "obs") df_obs_categorized else df_rpi_categorized) %>%
  filter(location == top_site)

# Order species by overall degree across all sites under this method
overall_poll_degree <- df_all %>%
  filter(method == top_method, ground_truth == 1) %>%
  group_by(higher_level) %>%
  summarise(overall_poll_degree = n_distinct(lower_level), .groups = "drop")

overall_plant_degree <- df_all %>%
  filter(method == top_method, ground_truth == 1) %>%
  group_by(lower_level) %>%
  summarise(overall_plant_degree = n_distinct(higher_level), .groups = "drop")

df_map <- df_site_cat %>%
  left_join(overall_poll_degree,  by = "higher_level") %>%
  left_join(overall_plant_degree, by = "lower_level") %>%
  mutate(
    link_display = recode(link_category,
      "locally_unique"   = "Locally unique",
      "unsupported"      = "Unsupported",
      "recurrent"        = "Recurrent",
      "model_elusive"    = "Model elusive",
      "likely_forbidden" = "Likely forbidden",
      "unconfirmed"      = "Unconfirmed",
      "locally_absent"   = "Locally absent",
      "possibly_missing" = "Possibly missing"
    )
  )

plant_order <- df_map %>%
  distinct(lower_level, overall_plant_degree) %>%
  arrange(desc(overall_plant_degree)) %>%
  pull(lower_level)

poll_order <- df_map %>%
  distinct(higher_level, overall_poll_degree) %>%
  arrange(desc(overall_poll_degree)) %>%
  pull(higher_level)

df_map <- df_map %>%
  mutate(
    lower_level  = factor(lower_level,  levels = plant_order),
    higher_level = factor(higher_level, levels = poll_order)
  )

map_interactive <- ggplot(df_map, aes(x = higher_level, y = lower_level,
                                      fill = link_display)) +
  geom_tile(color = "white") +
  scale_fill_brewer(palette = "Set2", name = "Link category") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 9, angle = 90, vjust = 0.5),
    axis.text.y  = element_text(size = 9),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 10)
  ) +
  labs(x = "Pollinator", y = "Plant") +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, " "), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  })) +
  scale_x_discrete(labels = function(x) lapply(strsplit(x, " "), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

map_interactive

## ---- 9. Map for most species-rich site ----
richness_table <- df_all %>%
  filter(ground_truth == 1) %>%
  group_by(location, method) %>%
  summarise(
    n_plants      = n_distinct(lower_level),
    n_pollinators = n_distinct(higher_level),
    n_species     = n_plants + n_pollinators,
    .groups = "drop"
  ) %>%
  arrange(desc(n_species))

print(richness_table)

rich_row    <- richness_table %>% slice_max(n_species, n = 1)
rich_site   <- rich_row$location
rich_method <- rich_row$method

cat(sprintf("\nMost species-rich site: %s (%s) — %d species\n",
            rich_site, rich_method, rich_row$n_species))

df_site_cat_rich <- (if (rich_method == "obs") df_obs_categorized else df_rpi_categorized) %>%
  filter(location == rich_site)

overall_poll_degree_rich <- df_all %>%
  filter(method == rich_method, ground_truth == 1) %>%
  group_by(higher_level) %>%
  summarise(overall_poll_degree = n_distinct(lower_level), .groups = "drop")

overall_plant_degree_rich <- df_all %>%
  filter(method == rich_method, ground_truth == 1) %>%
  group_by(lower_level) %>%
  summarise(overall_plant_degree = n_distinct(higher_level), .groups = "drop")

df_map_rich <- df_site_cat_rich %>%
  left_join(overall_poll_degree_rich,  by = "higher_level") %>%
  left_join(overall_plant_degree_rich, by = "lower_level") %>%
  mutate(
    link_display = recode(link_category,
      "locally_unique"   = "Locally unique",
      "unsupported"      = "Unsupported",
      "recurrent"        = "Recurrent",
      "model_elusive"    = "Model elusive",
      "likely_forbidden" = "Likely forbidden",
      "unconfirmed"      = "Unconfirmed",
      "locally_absent"   = "Locally absent",
      "possibly_missing" = "Possibly missing"
    )
  )

plant_order_rich <- df_map_rich %>%
  distinct(lower_level, overall_plant_degree) %>%
  arrange(desc(overall_plant_degree)) %>%
  pull(lower_level)

poll_order_rich <- df_map_rich %>%
  distinct(higher_level, overall_poll_degree) %>%
  arrange(desc(overall_poll_degree)) %>%
  pull(higher_level)

df_map_rich <- df_map_rich %>%
  mutate(
    lower_level  = factor(lower_level,  levels = plant_order_rich),
    higher_level = factor(higher_level, levels = poll_order_rich)
  )

map_rich <- ggplot(df_map_rich, aes(x = higher_level, y = lower_level,
                                    fill = link_display)) +
  geom_tile(color = "white") +
  scale_fill_brewer(palette = "Set2", name = "Link category") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 9, angle = 90, vjust = 0.5),
    axis.text.y  = element_text(size = 9),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 10)
  ) +
  labs(x = "Pollinator", y = "Plant") +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, " "), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  })) +
  scale_x_discrete(labels = function(x) lapply(strsplit(x, " "), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

map_rich
