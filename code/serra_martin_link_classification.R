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
col_confusion <- c("TP" = "lightsteelblue1",  "FP" = "lightsteelblue",
                   "TN" = "rosybrown",        "FN" = "rosybrown2")
col_validation <- c("Observed elsewhere"     = "sandybrown",
                    "Not observed elsewhere" = "thistle3")
col_subcats <- c(
  "recurrent"        = "coral",
  "locally_unique"   = "thistle3",
  "model_elusive"    = "coral2",
  "weak support"      = "thistle",
  "possibly_missing" = "coral1",
  "phantom"          = "thistle",
  "locally_absent"   = "coral3",
  "likely_forbidden" = "thistle2"
)
stratum_fill <- c(col_confusion, col_validation, col_subcats)

order_confusion  <- c("TP", "FP", "TN", "FN")
order_validation <- c("Not observed elsewhere", "Observed elsewhere")
order_subtypes   <- c(
  "locally_unique",   "phantom",
  "likely_forbidden", "weak support",
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
        is_unique     & original_binary == 1 & predicted_bin == 0 ~ "weak support",
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
  # "weak support" and "likely_forbidden" split by cross-method evidence;
  # "phantom" and all other subtypes pass through to L4 unchanged.
  order_L4 <- c(
    "locally_unique",
    "phantom",
    "weak support — have evidence", "weak support — no evidence",
    "recurrent", "possibly_missing", "locally_absent", "model_elusive",
    "likely_forbidden — have evidence", "likely_forbidden — no evidence"
  )
  col_L4 <- c(
    "locally_unique"                   = unname(col_subcats["locally_unique"]),
    "phantom"                          = unname(col_subcats["phantom"]),
    "weak support — have evidence"      = "plum2",    # weak support → more saturated
    "weak support — no evidence"        = "lavender", # weak support → lighter/greyer
    "recurrent"                        = unname(col_subcats["recurrent"]),
    "possibly_missing"                 = unname(col_subcats["possibly_missing"]),
    "locally_absent"                   = unname(col_subcats["locally_absent"]),
    "model_elusive"                    = unname(col_subcats["model_elusive"]),
    "likely_forbidden — have evidence" = "plum3",    # likely_forbidden → more saturated
    "likely_forbidden — no evidence"   = "thistle1"  # likely_forbidden → lighter/washed-out
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
        link_category == "weak support"      & add_obs  ~ "weak support — have evidence",
        link_category == "weak support"      & !add_obs ~ "weak support — no evidence",
        link_category == "likely_forbidden" & add_obs  ~ "likely_forbidden — have evidence",
        link_category == "likely_forbidden" & !add_obs ~ "likely_forbidden — no evidence",
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
# Level order mirrors the Sankey axis-4 grouping (used by both maps).
map_display_levels <- c(
  "Locally unique",
  "Recurrent",
  "Phantom — have evidence", "Phantom — no evidence",
  "Possibly missing",
  "Weak support",
  "Model elusive",
  "Likely forbidden — have evidence", "Likely forbidden — no evidence",
  "Locally absent"
)

# Palette grouped by confusion quadrant; hue signals biological meaning,
# lightness signals evidence.  Pairs sharing a hue (e.g. the two phantom or
# the two forbidden labels) use light = have evidence, dark = no evidence.
map_colors <- c(
  # TP (correctly predicted present) — teal
  "Locally unique"                    = "#76C7C0",  # light teal     (unique, TP)
  "Recurrent"                         = "#1F7E72",  # deep teal      (shared, TP)
  # FP (falsely predicted present) — warm coral/amber
  "Phantom — have evidence"           = "#F9C8B0",  # light coral    (all-zero FP, corroborated)
  "Phantom — no evidence"             = "#C05030",  # dark coral     (all-zero FP)
  "Possibly missing"                  = "#DBA040",  # amber          (obs-elsewhere FP)
  # FN (missed by model) — dusty rose
  "Weak support"                      = "#E8A8B8",  # light rose     (unique, FN)
  "Model elusive"                     = "#8A3050",  # deep rose      (shared, FN)
  # TN (correctly predicted absent) — lilac/purple
  "Likely forbidden — have evidence"  = "#D5C8F0",  # light lavender (all-zero TN, corroborated)
  "Likely forbidden — no evidence"    = "#5E3DA0",  # dark purple    (all-zero TN)
  "Locally absent"                    = "lightsteelblue"   # medium lilac   (obs-elsewhere TN)
)

df_site_cat <- (if (top_method == "obs") df_obs_categorized else df_rpi_categorized) %>%
  filter(location == top_site)

add_obs_ids_top <- if (top_method == "obs") add_obs_for_obs else add_obs_for_rpi

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
    add_obs = interaction_id %in% add_obs_ids_top,
    link_display = factor(case_when(
      link_category == "phantom"          & add_obs  ~ "Phantom — have evidence",
      link_category == "phantom"          & !add_obs ~ "Phantom — no evidence",
      link_category == "weak support"      ~ "Weak support",
      link_category == "likely_forbidden" & add_obs  ~ "Likely forbidden — have evidence",
      link_category == "likely_forbidden" & !add_obs ~ "Likely forbidden — no evidence",
      link_category == "locally_unique"   ~ "Locally unique",
      link_category == "recurrent"        ~ "Recurrent",
      link_category == "model_elusive"    ~ "Model elusive",
      link_category == "locally_absent"   ~ "Locally absent",
      link_category == "possibly_missing" ~ "Possibly missing"
    ), levels = map_display_levels)
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
  scale_fill_manual(values = map_colors, name = "Link category") +
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

add_obs_ids_rich <- if (rich_method == "obs") add_obs_for_obs else add_obs_for_rpi

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
    add_obs = interaction_id %in% add_obs_ids_rich,
    link_display = factor(case_when(
      link_category == "phantom"          & add_obs  ~ "Phantom — have evidence",
      link_category == "phantom"          & !add_obs ~ "Phantom — no evidence",
      link_category == "weak support"      ~ "Weak support",
      link_category == "likely_forbidden" & add_obs  ~ "Likely forbidden — have evidence",
      link_category == "likely_forbidden" & !add_obs ~ "Likely forbidden — no evidence",
      link_category == "locally_unique"   ~ "Locally unique",
      link_category == "recurrent"        ~ "Recurrent",
      link_category == "model_elusive"    ~ "Model elusive",
      link_category == "locally_absent"   ~ "Locally absent",
      link_category == "possibly_missing" ~ "Possibly missing"
    ), levels = map_display_levels)
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
  scale_fill_manual(values = map_colors, name = "Link category") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 13, angle = 90, vjust = 0.5),
    axis.text.y  = element_text(size = 13),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14)
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
# 
# ggsave("results/figures/richest_site_map_obs_validated_with_rpi.png",
#        plot = map_rich, width = 14, height = 7, dpi = 300, bg = "white")
# showtext::showtext_opts(dpi = 96)

# without axis titles
map_rich_blank <- ggplot(df_map_rich, aes(x = higher_level, y = lower_level,
                                    fill = link_display)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = map_colors, name = "Link category") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 13, angle = 90, vjust = 0.5),
    axis.text.y  = element_text(size = 13),
    axis.title = element_blank(),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14)
  ) +
  labs(x = "Pollinator", y = "Plant") +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, " "), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  })) +
  scale_x_discrete(labels = function(x) lapply(strsplit(x, " "), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

map_rich_blank

# ggsave("results/figures/map_rich_blank.svg",
#        plot = map_rich_blank, width = 14, height = 7, bg = "transparent")

## ---- 9b. Interactive Sankey (Plotly) ----------------------------------------
# Self-contained HTML — run the dedicated script from project root.
# Requires: plotly, htmlwidgets  (install.packages(c("plotly","htmlwidgets")))
source("interactive_sankey/make_sankey.R")

## ---- 9c. Map for richest rpi (camera) site, corroborated by observations ----
# Same logic as section 9 but pinned to method == "rpi" so the richest camera
# site is shown regardless of whether obs or rpi wins the global comparison.
# have evidence / no evidence: an interaction is corroborated when it also
# appears in the obs ground truth (add_obs_for_rpi).
rich_row_rpi   <- richness_table %>%
  dplyr::filter(method == "rpi") %>%
  dplyr::slice_max(n_species, n = 1)
rich_site_rpi  <- rich_row_rpi$location

cat(sprintf("\nRichest rpi site: %s — %d species\n",
            rich_site_rpi, rich_row_rpi$n_species))

df_site_cat_rpi_rich <- df_rpi_categorized %>%
  dplyr::filter(location == rich_site_rpi)

overall_poll_degree_rpi_rich <- df_all %>%
  dplyr::filter(method == "rpi", ground_truth == 1) %>%
  dplyr::group_by(higher_level) %>%
  dplyr::summarise(overall_poll_degree = n_distinct(lower_level), .groups = "drop")

overall_plant_degree_rpi_rich <- df_all %>%
  dplyr::filter(method == "rpi", ground_truth == 1) %>%
  dplyr::group_by(lower_level) %>%
  dplyr::summarise(overall_plant_degree = n_distinct(higher_level), .groups = "drop")

df_map_rpi_rich <- df_site_cat_rpi_rich %>%
  dplyr::left_join(overall_poll_degree_rpi_rich,  by = "higher_level") %>%
  dplyr::left_join(overall_plant_degree_rpi_rich, by = "lower_level") %>%
  dplyr::mutate(
    add_obs = interaction_id %in% add_obs_for_rpi,
    link_display = factor(dplyr::case_when(
      link_category == "phantom"          & add_obs  ~ "Phantom — have evidence",
      link_category == "phantom"          & !add_obs ~ "Phantom — no evidence",
      link_category == "weak support"      ~ "Weak support",
      link_category == "likely_forbidden" & add_obs  ~ "Likely forbidden — have evidence",
      link_category == "likely_forbidden" & !add_obs ~ "Likely forbidden — no evidence",
      link_category == "locally_unique"   ~ "Locally unique",
      link_category == "recurrent"        ~ "Recurrent",
      link_category == "model_elusive"    ~ "Model elusive",
      link_category == "locally_absent"   ~ "Locally absent",
      link_category == "possibly_missing" ~ "Possibly missing"
    ), levels = map_display_levels)
  )

plant_order_rpi_rich <- df_map_rpi_rich %>%
  dplyr::distinct(lower_level, overall_plant_degree) %>%
  dplyr::arrange(desc(overall_plant_degree)) %>%
  dplyr::pull(lower_level)

poll_order_rpi_rich <- df_map_rpi_rich %>%
  dplyr::distinct(higher_level, overall_poll_degree) %>%
  dplyr::arrange(desc(overall_poll_degree)) %>%
  dplyr::pull(higher_level)

df_map_rpi_rich <- df_map_rpi_rich %>%
  dplyr::mutate(
    lower_level  = factor(lower_level,  levels = plant_order_rpi_rich),
    higher_level = factor(higher_level, levels = poll_order_rpi_rich)
  )

map_rpi_rich <- ggplot(df_map_rpi_rich, aes(x = higher_level, y = lower_level,
                                             fill = link_display)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = map_colors, name = "Link category") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 11, angle = 90, vjust = 0.5),
    axis.text.y  = element_text(size = 11),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12)
  ) +
  labs(x = "Pollinator", y = "Plant",
       title = sprintf("Richest camera site: %s — rpi classifications, obs corroboration",
                       rich_site_rpi)) +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, " "), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  })) +
  scale_x_discrete(labels = function(x) lapply(strsplit(x, " "), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

map_rpi_rich

# ggsave("results/figures/map_richest_rpi_site.pdf",
#        plot = map_rpi_rich, width = 14, height = 7, device = cairo_pdf)
# ggsave("results/figures/map_richest_rpi_site.png",
#        plot = map_rpi_rich, width = 14, height = 7, dpi = 300, bg = "white")
## ---- 10. Sankey alluvial — ggsankey ----
# install once: devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)
# showtext + systemfonts: finds .otf fonts that extrafont misses.
library(showtext)
local({
  fam <- tryCatch(
    systemfonts::system_fonts() |>
      dplyr::filter(grepl("Amasis MT Pro", family, ignore.case = TRUE)),
    error = function(e) data.frame(path = character(0), style = character(0))
  )
  if (nrow(fam) == 0) {
    message("Amasis MT Pro not found. Available families containing 'Amasis': ",
            paste(unique(systemfonts::system_fonts()$family[grepl("amasis", systemfonts::system_fonts()$family, ignore.case = TRUE)]), collapse = ", "))
    return(invisible(NULL))
  }
  message("Amasis MT Pro found: ", paste(unique(fam$family), collapse = ", "))
  pick <- function(pat) {
    m <- fam |> dplyr::filter(grepl(pat, style, ignore.case = TRUE)) |> dplyr::slice(1)
    if (nrow(m)) m$path else fam$path[1]
  }
  sysfonts::font_add("AmasisMTPro",
                     regular = pick("Regular|Book|Roman|Light|Medium"),
                     bold    = pick("Bold"))
  showtext::showtext_auto()
})

make_sankey_validated <- function(df_categorized, add_obs_ids, method_label,
                                   validation_label = "Additional method",
                                   show_labels = TRUE,
                                   split_all = FALSE) {

  # ---- lavender gradient (locally_unique lightest → likely_forbidden darkest) -
  lav <- c(locally_unique    = "#C3BAD5",
            phantom           = "#AFA2C4",
            `weak support`    = "#9B8BB4",
            likely_forbidden  = "#8878A4")

  # ---- axis 4 ordering, colors, node IDs, and rank vectors -------------------
  # split_all = FALSE (default): only phantom, weak support, likely_forbidden
  #   are split into have / no evidence; all other categories pass through.
  # split_all = TRUE  (inspection): all 8 categories are split.
  if (!split_all) {

    order_L4 <- c(
      "locally_unique",
      "recurrent",
      "phantom — have evidence",         "phantom — no evidence",
      "possibly_missing",
      "weak support",
      "model_elusive",
      "likely_forbidden — have evidence","likely_forbidden — no evidence",
      "locally_absent"
    )

    col_L4 <- c(
      "locally_unique"                    = unname(lav["locally_unique"]),
      "recurrent"                         = unname(col_subcats["recurrent"]),
      "phantom — have evidence"           = "aquamarine",
      "phantom — no evidence"             = unname(lav["phantom"]),
      "possibly_missing"                  = unname(col_subcats["possibly_missing"]),
      "weak support"                      = unname(lav["weak support"]),
      "model_elusive"                     = unname(col_subcats["model_elusive"]),
      "likely_forbidden — have evidence"  = "aquamarine3",
      "likely_forbidden — no evidence"    = unname(lav["likely_forbidden"]),
      "locally_absent"                    = unname(col_subcats["locally_absent"])
    )

    node_id_vec <- c(
      TN = 1L, FN = 2L, FP = 3L, TP = 4L,
      "Not observed elsewhere" = 5L, "Observed elsewhere" = 6L,
      "likely_forbidden"                  =  7L,
      "likely_forbidden — no evidence"    =  8L,
      "likely_forbidden — have evidence"  =  9L,
      "weak support"                      = 10L,
      "phantom"                           = 13L,
      "phantom — no evidence"             = 14L,
      "phantom — have evidence"           = 15L,
      "locally_unique"                    = 16L,
      "locally_absent"                    = 17L,
      "model_elusive"                     = 18L,
      "possibly_missing"                  = 19L,
      "recurrent"                         = 20L
    )

    rank_vec <- c(
      "Not observed elsewhere" = 10, "Observed elsewhere" = 20,
      TN = 10, FN = 20, FP = 30, TP = 40,
      likely_forbidden = 10,
        "likely_forbidden — no evidence" = 11, "likely_forbidden — have evidence" = 12,
      "weak support" = 20,
      phantom = 30,
        "phantom — no evidence" = 31, "phantom — have evidence" = 32,
      locally_unique = 40, locally_absent = 50,
      model_elusive = 60, possibly_missing = 70, recurrent = 80
    )

    node_colors_l4 <- c(
      "phantom — have evidence"           = unname(col_L4["phantom — have evidence"]),
      "phantom — no evidence"             = unname(col_L4["phantom — no evidence"]),
      "likely_forbidden — have evidence"  = unname(col_L4["likely_forbidden — have evidence"]),
      "likely_forbidden — no evidence"    = unname(col_L4["likely_forbidden — no evidence"])
    )

  } else {
    # split_all = TRUE: all 8 categories split — aquamarine3 marks "have evidence"
    # consistently across all warm and cool category families.
    hev <- "aquamarine3"

    order_L4 <- c(
      "locally_unique — have evidence",    "locally_unique — no evidence",
      "recurrent — have evidence",         "recurrent — no evidence",
      "phantom — have evidence",           "phantom — no evidence",
      "possibly_missing — have evidence",  "possibly_missing — no evidence",
      "weak support — have evidence",      "weak support — no evidence",
      "model_elusive — have evidence",     "model_elusive — no evidence",
      "likely_forbidden — have evidence",  "likely_forbidden — no evidence",
      "locally_absent — have evidence",    "locally_absent — no evidence"
    )

    col_L4 <- c(
      "locally_unique — have evidence"    = hev,
      "locally_unique — no evidence"      = unname(lav["locally_unique"]),
      "recurrent — have evidence"         = hev,
      "recurrent — no evidence"           = unname(col_subcats["recurrent"]),
      "phantom — have evidence"           = hev,
      "phantom — no evidence"             = unname(lav["phantom"]),
      "possibly_missing — have evidence"  = hev,
      "possibly_missing — no evidence"    = unname(col_subcats["possibly_missing"]),
      "weak support — have evidence"      = hev,
      "weak support — no evidence"        = unname(lav["weak support"]),
      "model_elusive — have evidence"     = hev,
      "model_elusive — no evidence"       = unname(col_subcats["model_elusive"]),
      "likely_forbidden — have evidence"  = hev,
      "likely_forbidden — no evidence"    = unname(lav["likely_forbidden"]),
      "locally_absent — have evidence"    = hev,
      "locally_absent — no evidence"      = unname(col_subcats["locally_absent"])
    )

    # Axis-3 node IDs kept as-is (7–20); axis-4 split IDs start at 21.
    node_id_vec <- c(
      TN = 1L, FN = 2L, FP = 3L, TP = 4L,
      "Not observed elsewhere" = 5L, "Observed elsewhere" = 6L,
      "likely_forbidden"                  =  7L,
      "weak support"                      = 10L,
      "phantom"                           = 13L,
      "locally_unique"                    = 16L,
      "locally_absent"                    = 17L,
      "model_elusive"                     = 18L,
      "possibly_missing"                  = 19L,
      "recurrent"                         = 20L,
      "likely_forbidden — no evidence"    = 21L,
      "likely_forbidden — have evidence"  = 22L,
      "weak support — no evidence"        = 23L,
      "weak support — have evidence"      = 24L,
      "phantom — no evidence"             = 25L,
      "phantom — have evidence"           = 26L,
      "locally_unique — no evidence"      = 27L,
      "locally_unique — have evidence"    = 28L,
      "locally_absent — no evidence"      = 29L,
      "locally_absent — have evidence"    = 30L,
      "model_elusive — no evidence"       = 31L,
      "model_elusive — have evidence"     = 32L,
      "possibly_missing — no evidence"    = 33L,
      "possibly_missing — have evidence"  = 34L,
      "recurrent — no evidence"           = 35L,
      "recurrent — have evidence"         = 36L
    )

    rank_vec <- c(
      "Not observed elsewhere" = 10, "Observed elsewhere" = 20,
      TN = 10, FN = 20, FP = 30, TP = 40,
      likely_forbidden = 10,
        "likely_forbidden — no evidence" = 11, "likely_forbidden — have evidence" = 12,
      "weak support" = 20,
        "weak support — no evidence" = 21, "weak support — have evidence" = 22,
      phantom = 30,
        "phantom — no evidence" = 31, "phantom — have evidence" = 32,
      locally_unique = 40,
        "locally_unique — no evidence" = 41, "locally_unique — have evidence" = 42,
      locally_absent = 50,
        "locally_absent — no evidence" = 51, "locally_absent — have evidence" = 52,
      model_elusive = 60,
        "model_elusive — no evidence" = 61, "model_elusive — have evidence" = 62,
      possibly_missing = 70,
        "possibly_missing — no evidence" = 71, "possibly_missing — have evidence" = 72,
      recurrent = 80,
        "recurrent — no evidence" = 81, "recurrent — have evidence" = 82
    )

    node_colors_l4 <- c(
      "locally_unique — have evidence"    = hev,
      "locally_unique — no evidence"      = unname(lav["locally_unique"]),
      "recurrent — have evidence"         = hev,
      "recurrent — no evidence"           = unname(col_subcats["recurrent"]),
      "phantom — have evidence"           = hev,
      "phantom — no evidence"             = unname(lav["phantom"]),
      "possibly_missing — have evidence"  = hev,
      "possibly_missing — no evidence"    = unname(col_subcats["possibly_missing"]),
      "weak support — have evidence"      = hev,
      "weak support — no evidence"        = unname(lav["weak support"]),
      "model_elusive — have evidence"     = hev,
      "model_elusive — no evidence"       = unname(col_subcats["model_elusive"]),
      "likely_forbidden — have evidence"  = hev,
      "likely_forbidden — no evidence"    = unname(lav["likely_forbidden"]),
      "locally_absent — have evidence"    = hev,
      "locally_absent — no evidence"      = unname(col_subcats["locally_absent"])
    )
  }

  axis4_label <- paste0("Additional method: ", validation_label)
  axis_levels <- c("Within-network evaluation", "Contextual evidence",
                   "Link category", axis4_label)

  # ---- count links per unique path -------------------------------------------
  flows_tbl <- df_categorized %>%
    dplyr::filter(link_category != "unclassified", confusion != "UNK") %>%
    dplyr::mutate(
      add_obs = interaction_id %in% add_obs_ids,
      L4 = if (!split_all) {
        dplyr::case_when(
          link_category == "phantom"          & add_obs  ~ "phantom — have evidence",
          link_category == "phantom"          & !add_obs ~ "phantom — no evidence",
          link_category == "likely_forbidden" & add_obs  ~ "likely_forbidden — have evidence",
          link_category == "likely_forbidden" & !add_obs ~ "likely_forbidden — no evidence",
          TRUE ~ link_category
        )
      } else {
        dplyr::case_when(
          link_category == "locally_unique"   & add_obs  ~ "locally_unique — have evidence",
          link_category == "locally_unique"   & !add_obs ~ "locally_unique — no evidence",
          link_category == "recurrent"        & add_obs  ~ "recurrent — have evidence",
          link_category == "recurrent"        & !add_obs ~ "recurrent — no evidence",
          link_category == "phantom"          & add_obs  ~ "phantom — have evidence",
          link_category == "phantom"          & !add_obs ~ "phantom — no evidence",
          link_category == "possibly_missing" & add_obs  ~ "possibly_missing — have evidence",
          link_category == "possibly_missing" & !add_obs ~ "possibly_missing — no evidence",
          link_category == "weak support"     & add_obs  ~ "weak support — have evidence",
          link_category == "weak support"     & !add_obs ~ "weak support — no evidence",
          link_category == "model_elusive"    & add_obs  ~ "model_elusive — have evidence",
          link_category == "model_elusive"    & !add_obs ~ "model_elusive — no evidence",
          link_category == "likely_forbidden" & add_obs  ~ "likely_forbidden — have evidence",
          link_category == "likely_forbidden" & !add_obs ~ "likely_forbidden — no evidence",
          link_category == "locally_absent"   & add_obs  ~ "locally_absent — have evidence",
          link_category == "locally_absent"   & !add_obs ~ "locally_absent — no evidence",
          TRUE ~ link_category
        )
      }
    ) %>%
    dplyr::count(confusion, validation, link_category, L4, name = "value") %>%
    dplyr::mutate(L1 = confusion, L2 = validation, L3 = link_category, flow_color = L4)

  total_n <- sum(flows_tbl$value)

  # ---- reshape to long format ------------------------------------------------
  df_sankey <- flows_tbl %>%
    dplyr::select(L1, L2, L3, L4, value) %>%
    dplyr::mutate(alluvium = dplyr::row_number()) %>%
    tidyr::pivot_longer(c(L1, L2, L3, L4), names_to = "x_raw", values_to = "node") %>%
    dplyr::group_by(alluvium) %>%
    dplyr::arrange(x_raw, .by_group = TRUE) %>%
    dplyr::mutate(
      next_x_raw = dplyr::lead(x_raw),
      next_node  = dplyr::lead(node)
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
      node_clean = dplyr::if_else(
        node %in% c("TP", "FP", "TN", "FN"),
        node,
        sub("^(.)", "\\U\\1", stringr::str_replace_all(node, "_", " "), perl = TRUE)
      )
    ) %>%
    dplyr::select(-alluvium, -x_raw, -next_x_raw) %>%
    dplyr::group_by(x, node, next_x, next_node, node_clean) %>%
    dplyr::summarise(value = sum(value), .groups = "drop") %>%
    dplyr::mutate(
      .rank  = rank_vec[node],
      .nrank = rank_vec[next_node]
    ) %>%
    dplyr::group_by(x) %>%
    dplyr::arrange(.rank, .nrank, .by_group = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.rank, -.nrank) %>%
    dplyr::mutate(
      node_id      = node_id_vec[node],
      next_node_id = node_id_vec[next_node]
    )

  # ---- color palette ---------------------------------------------------------
  node_colors <- c(
    "TP"                     = unname(col_confusion["TP"]),
    "FP"                     = unname(col_confusion["FP"]),
    "FN"                     = unname(col_confusion["FN"]),
    "TN"                     = unname(col_confusion["TN"]),
    "Observed elsewhere"     = unname(col_validation["Observed elsewhere"]),
    "Not observed elsewhere" = unname(col_validation["Not observed elsewhere"]),
    "locally_unique"         = unname(lav["locally_unique"]),
    "recurrent"              = unname(col_subcats["recurrent"]),
    "phantom"                = unname(lav["phantom"]),
    "possibly_missing"       = unname(col_subcats["possibly_missing"]),
    "weak support"           = unname(lav["weak support"]),
    "model_elusive"          = unname(col_subcats["model_elusive"]),
    "likely_forbidden"       = unname(lav["likely_forbidden"]),
    "locally_absent"         = unname(col_subcats["locally_absent"]),
    node_colors_l4
  )

  # ---- node totals and label positions ---------------------------------------
  node_totals <- df_sankey %>%
    dplyr::group_by(x, node, node_clean) %>%
    dplyr::summarise(n_node = sum(value), .groups = "drop") %>%
    dplyr::mutate(pct_node = round(100 * n_node / total_n))

  p_probe <- ggplot(df_sankey,
                    aes(x = x, next_x = next_x, node = node_id,
                        next_node = next_node_id, value = value, fill = node)) +
    geom_sankey(space = 18) +
    geom_sankey_label(aes(label = node_id), space = 18)

  probe_layer <- ggplot_build(p_probe)$data[[2]]

  probe_y <- probe_layer %>%
    dplyr::transmute(
      x_int   = as.integer(round(n_x)),
      node_id = as.integer(node),
      y       = (ymin + ymax) / 2
    )

  bar_half     <- 0.002
  gap_side     <- 0.013
  mask_tile_df <- probe_layer %>%
    dplyr::transmute(
      x      = axis_levels[as.integer(round(n_x))],
      y_mid  = (ymin + ymax) / 2,
      height = ymax - ymin,
      width  = 2 * (bar_half + gap_side)
    )

  label_df <- node_totals %>%
    dplyr::mutate(
      node_id    = node_id_vec[node],
      x_int      = as.integer(x),
      node_label = paste0('bold("', node_clean, '")~plain("',
                          n_node, ' (', pct_node, '%)")')
    ) %>%
    dplyr::left_join(probe_y, by = c("x_int", "node_id")) %>%
    dplyr::arrange(x, node_id) %>%
    # axis 4 pass-through nodes (same name as axis 3) need no separate label;
    # in split_all mode every axis-4 node contains "—" so all are labelled.
    dplyr::filter(x_int < 4L | grepl("—", node))

  # ---- build the plot --------------------------------------------------------
  ggplot(df_sankey,
         aes(x = x, next_x = next_x, node = node_id, next_node = next_node_id,
             value = value, fill = node)) +
    geom_sankey(flow.alpha = flow_alpha, node.color = NA,
                node.width = 0.004, space = 18, smooth = 8) +
    geom_tile(
      data        = mask_tile_df,
      mapping     = aes(x = x, y = y_mid, width = width, height = height),
      fill        = "white",
      color       = NA,
      inherit.aes = FALSE
    ) +
    geom_sankey(flow.alpha = 0, node.color = NA,
                node.width = 0.004, space = 18, smooth = 8) +
    # ---- LABELS: set show_labels = FALSE in the function call to hide --------
    (if (show_labels) geom_text(
      data        = label_df,
      mapping     = aes(x = x, y = y, label = node_label),
      parse       = TRUE,
      nudge_x     = 0.08,
      hjust       = 0,
      size        = 3,
      color       = "grey15",
      inherit.aes = FALSE
    ) else NULL) +
    scale_fill_manual(values = node_colors, na.value = "grey80", guide = "none") +
    # ---- X EXPANSION: reduce right margin when labels are hidden -------------
    scale_x_discrete(expand = expansion(add = if (show_labels) c(0.3, 2.0) else c(0.3, 0.3))) +
    # ---- TITLES: set show_labels = FALSE to suppress title and axis text ----
    labs(title = if (show_labels) paste("Link classification — Sankey —", method_label) else NULL,
         subtitle = if (show_labels) paste("Within-network evaluation → Contextual evidence →",
                                           "Link category →", axis4_label) else NULL,
         x = NULL, y = NULL) +
    theme_sankey(base_size = 12, base_family = "AmasisMTPro") +
    theme(
      axis.text.x   = if (show_labels) element_text(size = 11, face = "bold") else element_blank(),
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

# ggsave("results/figures/sankey_obs_validated_with_cameras.pdf",
#        plot = gg_obs_sankey, width = 14, height = 6, device = cairo_pdf)
# showtext::showtext_opts(dpi = 300)
# ggsave("results/figures/sankey_obs_validated_with_cameras.png",
#        plot = gg_obs_sankey, width = 14, height = 6, dpi = 300, bg = "white")
# showtext::showtext_opts(dpi = 96)  # reset to screen dpi for viewer

# ---- unlabelled version (show_labels = FALSE) --------------------------------
gg_obs_sankey_clean <- make_sankey_validated(
  df_obs_categorized,
  add_obs_ids      = add_obs_for_obs,
  method_label     = "Direct observation (obs)",
  validation_label = "Raspberry Pi camera (rpi)",
  show_labels      = FALSE             # <-- toggle here
)
gg_obs_sankey_clean
# showtext::showtext_opts(dpi = 300)
# ggsave("results/figures/sankey_obs_validated_with_cameras_clean.pdf",
#        plot = gg_obs_sankey_clean, width = 14, height = 6, device = cairo_pdf)
# ggsave("results/figures/sankey_obs_validated_with_cameras_clean.png",
#        plot = gg_obs_sankey_clean, width = 14, height = 6, dpi = 300, bg = "white")
# showtext::showtext_opts(dpi = 96)

gg_rpi_sankey <- make_sankey_validated(
  df_rpi_categorized,
  add_obs_ids      = add_obs_for_rpi,
  method_label     = "Raspberry Pi camera (rpi)",
  validation_label = "Direct observation (obs)"
)
gg_rpi_sankey
# ggsave("results/figures/sankey_rpi_validated_with_obs.pdf",
#        plot = gg_rpi_sankey, width = 14, height = 6, device = cairo_pdf)
# showtext::showtext_opts(dpi = 300)
# ggsave("results/figures/sankey_rpi_validated_with_obs.png",
#        plot = gg_rpi_sankey, width = 14, height = 6, dpi = 300, bg = "white")
# showtext::showtext_opts(dpi = 96)

# ---- unlabelled version (show_labels = FALSE) --------------------------------
gg_rpi_sankey_clean <- make_sankey_validated(
  df_rpi_categorized,
  add_obs_ids      = add_obs_for_rpi,
  method_label     = "Raspberry Pi camera (rpi)",
  validation_label = "Direct observation (obs)",
  show_labels      = FALSE             # <-- toggle here
)
gg_rpi_sankey_clean
# showtext::showtext_opts(dpi = 300)
# ggsave("results/figures/sankey_rpi_validated_with_obs_clean.pdf",
#        plot = gg_rpi_sankey_clean, width = 14, height = 6, device = cairo_pdf)
# ggsave("results/figures/sankey_rpi_validated_with_obs_clean.png",
#        plot = gg_rpi_sankey_clean, width = 14, height = 6, dpi = 300, bg = "white")
# showtext::showtext_opts(dpi = 96)

# ---- inspection versions: all 8 categories split into have / no evidence ----
# split_all = TRUE adds have evidence / no evidence at axis 4 for every category,
# not just phantom / weak support / likely_forbidden.  For inspection only.
gg_obs_sankey_full <- make_sankey_validated(
  df_obs_categorized,
  add_obs_ids      = add_obs_for_obs,
  method_label     = "Direct observation (obs) — all categories",
  validation_label = "Raspberry Pi camera (rpi)",
  split_all        = TRUE
)
gg_obs_sankey_full

gg_rpi_sankey_full <- make_sankey_validated(
  df_rpi_categorized,
  add_obs_ids      = add_obs_for_rpi,
  method_label     = "Raspberry Pi camera (rpi) — all categories",
  validation_label = "Direct observation (obs)",
  split_all        = TRUE
)
gg_rpi_sankey_full

## ---- 11. Summary: unobserved link breakdown ----
# For a given sampling method, pooling all six sites:
#
# LINK-SITE RECORDS: the unit of observation.  Each record = one unique species-pair
#   (interaction) assessed at one site under this method.  The same species pair
#   can appear in up to six records (once per site), so the total number of
#   link-site records is always larger than the number of unique species pairs.
#
# UNIQUE SPECIES INTERACTIONS: distinct species-pair combinations in the dataset,
#   regardless of how many sites they were assessed at.  A pair that was assessed
#   at all six sites contributes six link-site records but only one unique interaction.
#
# UNOBSERVED RECORDS are split by the model's prediction:
#   FP (false positive)  — unobserved at this site, yet the model predicted present.
#     · phantom          — never recorded at ANY site under this method
#                          (the interaction may be genuinely impossible, or simply missed)
#     · possibly_missing — recorded at ≥1 other site under this method,
#                          but absent here (likely a detection gap)
#   TN (true negative)   — unobserved at this site and the model correctly predicted absent.
#     · likely_forbidden — never recorded at ANY site under this method
#                          (the model concurs: the interaction is probably impossible)
#     · locally_absent   — recorded at ≥1 other site under this method,
#                          but model correctly predicts absence here (locally unsuitable)
#
# CROSS-METHOD EVIDENCE ("have evidence"): the interaction_id appears as
#   ground_truth == 1 in the orthogonal sampling method (cameras for obs, or
#   vice versa).  "% corroborated" = have evidence / category total.
summarise_unobserved <- function(df_cat, add_obs_ids, method_name, other_method_name) {
  df <- df_cat %>%
    dplyr::filter(link_category != "unclassified", confusion != "UNK") %>%
    dplyr::mutate(add_obs = interaction_id %in% add_obs_ids)

  total   <- nrow(df)
  n_obs   <- sum(df$original_binary == 1)
  n_unobs <- total - n_obs
  n_fp    <- sum(df$confusion == "FP")
  n_tn    <- sum(df$confusion == "TN")

  n_ids_total <- dplyr::n_distinct(df$interaction_id)
  n_ids_never <- df %>%
    dplyr::filter(is_all_zero) %>%
    dplyr::pull(interaction_id) %>%
    dplyr::n_distinct()
  n_ids_seen  <- n_ids_total - n_ids_never

  make_tbl <- function(confusion_label, n_denom) {
    df %>%
      dplyr::filter(confusion == confusion_label) %>%
      dplyr::group_by(link_category) %>%
      dplyr::summarise(
        have_evidence = sum(add_obs),
        no_evidence   = sum(!add_obs),
        total_cat     = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        pct_of_group = sprintf("%.1f%%", 100 * total_cat  / n_denom),
        pct_hev      = sprintf("%.1f%%", 100 * have_evidence / total_cat)
      ) %>%
      dplyr::arrange(dplyr::desc(total_cat))
  }

  fp_tbl <- make_tbl("FP", n_fp)
  tn_tbl <- make_tbl("TN", n_tn)

  print_tbl <- function(tbl, denom, denom_label) {
    hdr <- sprintf("  %-20s  %7s  %-14s  %-16s  %s\n",
                   "Category", "n",
                   sprintf("%% of %s", denom_label),
                   sprintf("have evidence (%s)", other_method_name),
                   "% corroborated")
    cat(hdr)
    cat(strrep("-", nchar(hdr) - 1L), "\n")
    for (i in seq_len(nrow(tbl))) {
      r <- tbl[i, ]
      cat(sprintf("  %-20s  %7d  %-14s  %-16d  %s\n",
                  r$link_category, r$total_cat, r$pct_of_group,
                  r$have_evidence, r$pct_hev))
    }
    cat(strrep("-", nchar(hdr) - 1L), "\n")
    cat(sprintf("  %-20s  %7d  %-14s  %d\n",
                "Total", denom, "100.0%", sum(tbl$have_evidence)))
  }

  sep <- strrep("=", 72)
  cat(sprintf("\n%s\n", sep))
  cat(sprintf("  METHOD: %s\n", method_name))
  cat(sprintf("%s\n", sep))

  cat(sprintf("\n--- LINK-SITE RECORDS (each species pair x site = 1 record) ---\n\n"))
  cat(sprintf("  Total link-site records:                          %d\n", total))
  cat(sprintf("  Observed (ground truth = 1):                      %d  (%.1f%% of %d total records)\n",
              n_obs,   100 * n_obs   / total, total))
  cat(sprintf("  Unobserved (ground truth = 0):                    %d  (%.1f%% of %d total records)\n",
              n_unobs, 100 * n_unobs / total, total))
  cat(sprintf("    of which predicted present — FP:                %d  (%.1f%% of %d unobserved records)\n",
              n_fp, 100 * n_fp / n_unobs, n_unobs))
  cat(sprintf("    of which predicted absent  — TN:                %d  (%.1f%% of %d unobserved records)\n",
              n_tn, 100 * n_tn / n_unobs, n_unobs))

  cat(sprintf("\n--- UNIQUE SPECIES INTERACTIONS (pooled across all sites) ---\n\n"))
  cat(sprintf("  Total unique species pairs in dataset:            %d\n", n_ids_total))
  cat(sprintf("  Observed at ≥1 site (ground truth = 1 anywhere): %d  (%.1f%% of %d unique pairs)\n",
              n_ids_seen,  100 * n_ids_seen  / n_ids_total, n_ids_total))
  cat(sprintf("  Never observed at any site:                       %d  (%.1f%% of %d unique pairs)\n",
              n_ids_never, 100 * n_ids_never / n_ids_total, n_ids_total))

  cat(sprintf(
    "\n--- UNOBSERVED & PREDICTED PRESENT (FP, n = %d records) ---\n", n_fp))
  cat(sprintf(
    "    phantom        = predicted present, never recorded at any site under %s\n",
    method_name))
  cat(sprintf(
    "    possibly_missing = predicted present, recorded at ≥1 other site under %s\n\n",
    method_name))
  print_tbl(fp_tbl, n_fp, "FP")

  cat(sprintf(
    "\n--- UNOBSERVED & PREDICTED ABSENT (TN, n = %d records) ---\n", n_tn))
  cat(sprintf(
    "    likely_forbidden = predicted absent, never recorded at any site under %s\n",
    method_name))
  cat(sprintf(
    "    locally_absent   = predicted absent, recorded at ≥1 other site under %s\n\n",
    method_name))
  print_tbl(tn_tbl, n_tn, "TN")
  cat("\n")

  invisible(list(fp = fp_tbl, tn = tn_tbl))
}

summarise_unobserved(
  df_obs_categorized,
  add_obs_ids      = add_obs_for_obs,
  method_name      = "Direct observation (obs)",
  other_method_name = "rpi"
)

summarise_unobserved(
  df_rpi_categorized,
  add_obs_ids      = add_obs_for_rpi,
  method_name      = "Raspberry Pi camera (rpi)",
  other_method_name = "obs"
)
