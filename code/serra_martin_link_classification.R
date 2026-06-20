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
  "low support"      = "thistle",
  "possibly_missing" = "coral2",
  "phantom"          = "thistle",
  "locally_absent"   = "coral1",
  "likely_forbidden" = "thistle2"
)
stratum_fill <- c(col_confusion, col_validation, col_subcats)

order_confusion  <- c("TP", "FP", "TN", "FN")
order_validation <- c("Not observed elsewhere", "Observed elsewhere")
order_subtypes   <- c(
  "locally_unique",   "phantom",
  "likely_forbidden", "low support",
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
        is_unique     & original_binary == 1 & predicted_bin == 0 ~ "low support",
        is_shared     & original_binary == 1 & predicted_bin == 1 ~ "recurrent",
        is_shared     & original_binary == 1 & predicted_bin == 0 ~ "model_elusive",
        is_all_zero   & original_binary == 0 & predicted_bin == 0 ~ "likely_forbidden",
        is_all_zero   & original_binary == 0 & predicted_bin == 1 ~ "phantom",
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
  # "low support" and "likely_forbidden" split by cross-method evidence;
  # "phantom" and all other subtypes pass through to L4 unchanged.
  order_L4 <- c(
    "locally_unique",
    "phantom",
    "low support — Have evidence", "low support — No evidence",
    "recurrent", "possibly_missing", "locally_absent", "model_elusive",
    "likely_forbidden — Have evidence", "likely_forbidden — No evidence"
  )
  col_L4 <- c(
    "locally_unique"                   = unname(col_subcats["locally_unique"]),
    "phantom"                          = unname(col_subcats["phantom"]),
    "low support — Have evidence"      = "plum2",    # low support → more saturated
    "low support — No evidence"        = "lavender", # low support → lighter/greyer
    "recurrent"                        = unname(col_subcats["recurrent"]),
    "possibly_missing"                 = unname(col_subcats["possibly_missing"]),
    "locally_absent"                   = unname(col_subcats["locally_absent"]),
    "model_elusive"                    = unname(col_subcats["model_elusive"]),
    "likely_forbidden — Have evidence" = "plum3",    # likely_forbidden → more saturated
    "likely_forbidden — No evidence"   = "thistle1"  # likely_forbidden → lighter/washed-out
  )
  stratum_fill_4    <- c(stratum_fill, col_L4)
  axis4_label       <- validation_label
  order_axis4       <- c("Confusion", "Validation", "Subtype", axis4_label)
  label_color_map_4 <- c(label_color_map, setNames("mistyrose4", axis4_label))

  flows_tbl <- df_categorized %>%
    filter(link_category != "unclassified", confusion != "UNK") %>%
    mutate(
      add_obs = interaction_id %in% add_obs_ids,
      L4 = dplyr::case_when(
        link_category == "low support"      & add_obs  ~ "low support — Have evidence",
        link_category == "low support"      & !add_obs ~ "low support — No evidence",
        link_category == "likely_forbidden" & add_obs  ~ "likely_forbidden — Have evidence",
        link_category == "likely_forbidden" & !add_obs ~ "likely_forbidden — No evidence",
        TRUE ~ link_category
      )
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
    ) %>%
    dplyr::group_by(alluvium_id) %>%
    dplyr::arrange(axis, .by_group = TRUE) %>%
    dplyr::mutate(
      dest_fill = coalesce(as.character(dplyr::lead(stratum)), as.character(stratum))
    ) %>%
    dplyr::ungroup()

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
  add_obs_ids      = add_obs_for_obs, # these are validations using the orthogonal method (unique interactions that are ground_truth == 1 using cameras, in this case)
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
      "low support"      = "Low support",
      "recurrent"        = "Recurrent",
      "model_elusive"    = "Model elusive",
      "likely_forbidden" = "Likely forbidden",
      "phantom"      = "Phantom",
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
      "low support"      = "Low support",
      "recurrent"        = "Recurrent",
      "model_elusive"    = "Model elusive",
      "likely_forbidden" = "Likely forbidden",
      "phantom"      = "Phantom",
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

# pdf(file   = "results/figures/map_richest_site.pdf",
#     width  = 14,    # inches
#     height = 7,
#     family = "Helvetica"   # or another installed font
# )
# map_rich
# dev.off()

## ---- 10. Sankey alluvial — ggsankey ----
# install once: devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)

make_sankey_validated <- function(df_categorized, add_obs_ids, method_label,
                                   validation_label = "Additional method") {

  # ---- axis 4 ordering and colors --------------------------------------------
  # The fourth axis shows 10 final categories: low support and likely_forbidden
  # are each split into "Have evidence" (seen by the other method) and
  # "No evidence" (not corroborated). All other categories pass through unchanged.
  # L4 categories grouped by confusion pair (TP→, FP→, FN→, TN→) so that
  # axis 3 → axis 4 ribbons stay within each confusion group and do not cross.
  order_L4 <- c(
    "locally_unique",
    "recurrent",
    "phantom — Have evidence", "phantom — No evidence",
    "possibly_missing",
    "low support — Have evidence", "low support — No evidence",
    "model_elusive",
    "likely_forbidden — Have evidence", "likely_forbidden — No evidence",
    "locally_absent"
  )

  # "Have evidence" variants: saturated plum tones.
  # "No evidence" variants: medium thistle tones, consistent with the FN/TN family.
  col_L4 <- c(
    "locally_unique"                   = unname(col_subcats["locally_unique"]),
    "recurrent"                        = unname(col_subcats["recurrent"]),
    "phantom — Have evidence"          = "plum1",
    "phantom — No evidence"            = "thistle",
    "possibly_missing"                 = unname(col_subcats["possibly_missing"]),
    "low support — Have evidence"      = "plum2",
    "low support — No evidence"        = "thistle2",
    "model_elusive"                    = unname(col_subcats["model_elusive"]),
    "likely_forbidden — Have evidence" = "plum3",
    "likely_forbidden — No evidence"   = "thistle3",
    "locally_absent"                   = unname(col_subcats["locally_absent"])
  )

  axis4_label <- paste0("Additional method: ", validation_label)
  axis_levels <- c("Within-network evaluation", "Contextual evidence",
                   "Link category", axis4_label)

  # ---- count links per unique path -------------------------------------------
  # flows_tbl has one row per combination of (confusion × validation × subtype × L4).
  # value = number of individual links (summed across all 6 sites) on that path.
  # This is what you see in the diagram: e.g. TP = all correctly-predicted present
  # links across every site for this sampling method.
  flows_tbl <- df_categorized %>%
    dplyr::filter(link_category != "unclassified", confusion != "UNK") %>%
    dplyr::mutate(
      add_obs = interaction_id %in% add_obs_ids,
      # assign final L4 category based on cross-method evidence
      L4 = dplyr::case_when(
        link_category == "phantom"          & add_obs  ~ "phantom — Have evidence",
        link_category == "phantom"          & !add_obs ~ "phantom — No evidence",
        link_category == "low support"      & add_obs  ~ "low support — Have evidence",
        link_category == "low support"      & !add_obs ~ "low support — No evidence",
        link_category == "likely_forbidden" & add_obs  ~ "likely_forbidden — Have evidence",
        link_category == "likely_forbidden" & !add_obs ~ "likely_forbidden — No evidence",
        TRUE ~ link_category
      )
    ) %>%
    dplyr::count(confusion, validation, link_category, L4, name = "value") %>%
    dplyr::mutate(L1 = confusion, L2 = validation, L3 = link_category, flow_color = L4)

  # total links across all categories — denominator for percentage labels
  total_n <- sum(flows_tbl$value)

  # ---- reshape to long format (one row per axis per path) --------------------
  # pivot_longer turns the 10-path wide table into 40 rows (10 × 4 axes).
  # lead() gives each row a pointer to the next axis's node so ggsankey can
  # draw the ribbon connecting them.
  # Integer stacking-order IDs: 1 = bottom of the diagram, highest = top.
  # StatSankeyFlow sorts by the node aesthetic value; using integers ensures
  # numerical (not alphabetical) ordering so bars stack in our desired order.
  node_id_vec <- c(
    TN = 1L, FN = 2L, FP = 3L, TP = 4L,
    "Not observed elsewhere" = 5L, "Observed elsewhere" = 6L,
    "likely_forbidden"                 =  7L,
    "likely_forbidden — No evidence"   =  8L,
    "likely_forbidden — Have evidence" =  9L,
    "low support"                      = 10L,
    "low support — No evidence"        = 11L,
    "low support — Have evidence"      = 12L,
    "phantom"                          = 13L,
    "phantom — No evidence"            = 14L,
    "phantom — Have evidence"          = 15L,
    "locally_unique"                   = 16L,
    "locally_absent"                   = 17L,
    "model_elusive"                    = 18L,
    "possibly_missing"                 = 19L,
    "recurrent"                        = 20L
  )

  df_sankey <- flows_tbl %>%
    dplyr::select(L1, L2, L3, L4, value) %>%
    dplyr::mutate(alluvium = dplyr::row_number()) %>%
    tidyr::pivot_longer(c(L1, L2, L3, L4), names_to = "x_raw", values_to = "node") %>%
    dplyr::group_by(alluvium) %>%
    dplyr::arrange(x_raw, .by_group = TRUE) %>%
    dplyr::mutate(
      next_x_raw = dplyr::lead(x_raw),  # which axis comes next
      next_node  = dplyr::lead(node)    # what category is at that next axis
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      x      = dplyr::recode(x_raw,
                 L1 = "Within-network evaluation", L2 = "Contextual evidence",
                 L3 = "Link category",             L4 = axis4_label),
      next_x = dplyr::recode(next_x_raw,
                 L1 = "Within-network evaluation", L2 = "Contextual evidence",
                 L3 = "Link category",             L4 = axis4_label),
      x      = factor(x,      levels = axis_levels),
      next_x = factor(next_x, levels = axis_levels),
      # TP/FP/TN/FN stay uppercase; all other names get sentence case
      node_clean = dplyr::if_else(
        node %in% c("TP", "FP", "TN", "FN"),
        node,
        sub("^(.)", "\\U\\1", stringr::str_replace_all(node, "_", " "), perl = TRUE)
      )
    ) %>%
    dplyr::select(-alluvium, -x_raw, -next_x_raw) %>%
    # collapse rows that share the same (source → destination) edge into one ribbon
    dplyr::group_by(x, node, next_x, next_node, node_clean) %>%
    dplyr::summarise(value = sum(value), .groups = "drop") %>%
    # ---- control stacking order -------------------------------------------------
    # ggsankey sorts nodes by the node aesthetic value (alphabetical for
    # characters, ascending for numbers).  We add a numeric node_id column
    # (1 = bottom, 20 = top) to force exact ordering — see node_id_vec below.
    # These .rank / .nrank weights are NOT link counts; they are just sorting
    # keys used to arrange rows consistently before that step.
    # Desired top-to-bottom layout:
    #   Axis 1: TP, FP, FN, TN
    #   Axis 2: Observed elsewhere, Not observed elsewhere
    #   Axis 3: recurrent / possibly_missing / model_elusive / locally_absent
    #           (Observed group), then locally_unique / phantom / low support /
    #           likely_forbidden (Not observed group)
    dplyr::mutate(
      .rank = c(
        "Not observed elsewhere" = 10, "Observed elsewhere" = 20,
        TN = 10, FN = 20, FP = 30, TP = 40,
        likely_forbidden = 10, "likely_forbidden — No evidence" = 11,
                               "likely_forbidden — Have evidence" = 12,
        "low support" = 20, "low support — No evidence" = 21,
                            "low support — Have evidence" = 22,
        phantom = 30, "phantom — No evidence" = 31,
                      "phantom — Have evidence" = 32,
        locally_unique = 40,
        locally_absent = 50,
        model_elusive  = 60,
        possibly_missing = 70,
        recurrent = 80
      )[node],
      .nrank = c(
        "Not observed elsewhere" = 10, "Observed elsewhere" = 20,
        TN = 10, FN = 20, FP = 30, TP = 40,
        likely_forbidden = 10, "likely_forbidden — No evidence" = 11,
                               "likely_forbidden — Have evidence" = 12,
        "low support" = 20, "low support — No evidence" = 21,
                            "low support — Have evidence" = 22,
        phantom = 30, "phantom — No evidence" = 31,
                      "phantom — Have evidence" = 32,
        locally_unique = 40,
        locally_absent = 50,
        model_elusive  = 60,
        possibly_missing = 70,
        recurrent = 80
      )[next_node]
    ) %>%
    dplyr::group_by(x) %>%
    dplyr::arrange(.rank, .nrank, .by_group = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.rank, -.nrank) %>%
    # Map each node name to its integer stacking ID (1 = bottom, 20 = top).
    # The ggsankey stats sort by these integers numerically, giving the
    # correct bottom-to-top bar/label order on every axis.
    dplyr::mutate(
      node_id      = node_id_vec[node],
      next_node_id = node_id_vec[next_node]
    )

  # ---- color palette: one entry per node name ---------------------------------
  # fill = node keeps every source-node's fill constant across its outgoing
  # ribbons, which prevents ggsankey from inserting extra spacing between
  # ribbons of the same node when they carry different fill values.
  # Edit any entry here to recolor that node (bar + outgoing ribbon) everywhere
  # it appears in the diagram.  Axis 4 split-nodes get their own entries.
  node_colors <- c(
    # Axis 1 — Within-network evaluation (confusion matrix)
    "TP"                               = unname(col_confusion["TP"]),
    "FP"                               = unname(col_confusion["FP"]),
    "FN"                               = unname(col_confusion["FN"]),
    "TN"                               = unname(col_confusion["TN"]),
    # Axis 2 — Contextual evidence
    "Observed elsewhere"               = unname(col_validation["Observed elsewhere"]),
    "Not observed elsewhere"           = unname(col_validation["Not observed elsewhere"]),
    # Axis 3 — Link category (also covers singleton Axis 4 nodes of the same name)
    "locally_unique"                   = unname(col_subcats["locally_unique"]),
    "recurrent"                        = unname(col_subcats["recurrent"]),
    "phantom"                          = unname(col_subcats["phantom"]),
    "possibly_missing"                 = unname(col_subcats["possibly_missing"]),
    "low support"                      = unname(col_subcats["low support"]),
    "model_elusive"                    = unname(col_subcats["model_elusive"]),
    "likely_forbidden"                 = "aquamarine3",
    "locally_absent"                   = unname(col_subcats["locally_absent"]),
    # Axis 4 — Additional method subcategories (split nodes only)
    "phantom — Have evidence"          = unname(col_L4["phantom — Have evidence"]),
    "phantom — No evidence"            = unname(col_L4["phantom — No evidence"]),
    "low support — Have evidence"      = unname(col_L4["low support — Have evidence"]),
    "low support — No evidence"        = unname(col_L4["low support — No evidence"]),
    "likely_forbidden — Have evidence" = "aquamarine4",
    "likely_forbidden — No evidence"   = unname(col_L4["likely_forbidden — No evidence"])
  )

  # ---- compute per-node counts and percentages for labels --------------------
  # For each (axis, category), sum all flow values passing through that node
  # to get how many links it represents and what share of the total that is.
  node_totals <- df_sankey %>%
    dplyr::group_by(x, node, node_clean) %>%
    dplyr::summarise(n_node = sum(value), .groups = "drop") %>%
    dplyr::mutate(pct_node = round(100 * n_node / total_n))

  # StatSankeyText aligns labels with bars correctly but sorts nodes
  # alphabetically — wrong order for our diagram.  Build a probe ggplot with
  # geom_sankey_label and extract the bar-midpoint y-coordinates from
  # ggplot_build(), keyed by node name.  Join these onto node_totals (sorted by
  # node_id) so geom_text gets ggsankey's own positions in our desired order.
  p_probe <- ggplot(df_sankey,
                    aes(x = x, next_x = next_x, node = node_id,
                        next_node = next_node_id, value = value, fill = node)) +
    geom_sankey(space = 8) +
    geom_sankey_label(aes(label = node_id), space = 8)

  # The ggsankey `node` aesthetic (= node_id integer) survives ggplot_build
  # unchanged and uniquely identifies each bar.  Use it as the join key so we
  # don't depend on label or fill column names, which vary across ggsankey versions.
  probe_y <- ggplot_build(p_probe)$data[[2]] %>%
    dplyr::transmute(
      x_int   = as.integer(round(n_x)),
      node_id = as.integer(node),
      y       = (ymin + ymax) / 2
    )

  label_df <- node_totals %>%
    dplyr::mutate(
      node_id    = node_id_vec[node],
      x_int      = as.integer(x),
      node_label = paste0('bold("', node_clean, '")~plain("',
                          n_node, ' (', pct_node, '%)")')
    ) %>%
    dplyr::left_join(probe_y, by = c("x_int", "node_id")) %>%
    dplyr::arrange(x, node_id)

  # ---- build the plot --------------------------------------------------------
  ggplot(df_sankey,
         aes(x = x, next_x = next_x, node = node_id, next_node = next_node_id,
             value = value, fill = node)) +
    geom_sankey(flow.alpha = flow_alpha, node.color = "white",
                node.width = 0.008, space = 8, smooth = 8) +
    geom_text(
      data        = label_df,
      mapping     = aes(x = x, y = y, label = node_label),
      parse       = TRUE,
      nudge_x     = 0.05,
      hjust       = 0,
      size        = 3,
      color       = "grey30",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = node_colors, na.value = "grey80", guide = "none") +
    scale_x_discrete(expand = expansion(add = c(0.3, 2.0))) +
    labs(
      title    = paste("Link classification — Sankey —", method_label),
      subtitle = paste("Within-network evaluation → Contextual evidence →",
                       "Link category →", axis4_label),
      x = NULL, y = NULL
    ) +
    theme_sankey(base_size = 12) +
    theme(
      axis.text.x   = element_text(size = 11, face = "bold"),
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30"),
      plot.margin   = margin(20, 20, 20, 20)
    )
}

gg_obs_sankey <- make_sankey_validated(
  df_obs_categorized,
  add_obs_ids      = add_obs_for_obs,
  method_label     = "Direct observation (obs)",
  validation_label = "Raspberry Pi camera (rpi)"
)
gg_obs_sankey

# pdf(file   = "results/figures/sankey_obs_validated_with_cameras.pdf",
#     width  = 12,    # inches
#     height = 6,
#     family = "Helvetica"   # or another installed font
# )
# gg_obs_sankey
# dev.off()

gg_rpi_sankey <- make_sankey_validated(
  df_rpi_categorized,
  add_obs_ids      = add_obs_for_rpi,
  method_label     = "Raspberry Pi camera (rpi)",
  validation_label = "Direct observation (obs)"
)
gg_rpi_sankey
