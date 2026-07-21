## ---- Cross-location cross-method link corroboration ----
# Corroborates eco_ILP link predictions using spatial and methodological evidence.
# Data sources:
#   Section 1 — Dryad plant-pollinator networks (Arstingstall et al. 2021, Mol. Ecol.)
#               3 locations × 2 methods (observation, DNA metabarcoding)
#   Section 2 — FrugInt frugivory networks (MN_2024, mist-netting)
#               2 Pistacia sites (Hato Ratón, Southern site) for spatial corroboration
# Predictions from eco_ILP (full matrix, no link withholding, threshold = 0.5).

# ---- load libraries ----
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggalluvial)
library(scales)
library(patchwork)

# ---- Section 1: Dryad plant-pollinator — load prediction data ----
Starkey_metabarcoding <- read_csv("results/predictions/Starkey_metabarcoding_prediction_results.csv", show_col_types = FALSE)
Starkey_observations <- read_csv("results/predictions/Starkey_observation_prediction_results.csv", show_col_types = FALSE)
Threemile_metabarcoding <- read_csv("results/predictions/Threemile_metabarcoding_prediction_results.csv", show_col_types = FALSE)
Threemile_observations <- read_csv("results/predictions/Threemile_observation_prediction_results.csv", show_col_types = FALSE)
Zumwalt_metabarcoding <- read_csv("results/predictions/Zumwalt_metabarcoding_prediction_results.csv", show_col_types = FALSE)
Zumwalt_observations <- read_csv("results/predictions/Zumwalt_observation_prediction_results.csv", show_col_types = FALSE)
metaweb_metabarcoding <- read_csv("results/predictions/metaweb_metabarcoding_prediction_results.csv", show_col_types = FALSE)
metaweb_observations <- read_csv("results/predictions/metaweb_observation_prediction_results.csv", show_col_types = FALSE)


Starkey_metabarcoding <- Starkey_metabarcoding %>% mutate(location = "Starkey",
                                                          method = "metabarcoding")
Starkey_observations <- Starkey_observations %>% mutate(location = "Starkey",
                                                          method = "observations")
Threemile_metabarcoding <- Threemile_metabarcoding %>% mutate(location = "Threemile",
                                                          method = "metabarcoding")
Threemile_observations <- Threemile_observations %>% mutate(location = "Threemile",
                                                        method = "observations")
Zumwalt_metabarcoding <- Zumwalt_metabarcoding %>% mutate(location = "Zumwalt",
                                                              method = "metabarcoding")
Zumwalt_observations <- Zumwalt_observations %>% mutate(location = "Zumwalt",
                                                            method = "observations")

df_all <- bind_rows(Starkey_metabarcoding, Starkey_observations, Threemile_metabarcoding,
                    Threemile_observations, Zumwalt_metabarcoding, Zumwalt_observations) %>%
  mutate(
    interaction_id = paste(higher_level, lower_level, sep = "___")
  )

# ---- diagnostic checks ----
# These checks verify that interaction IDs match across locations and that
# the "observed elsewhere" logic can actually find shared links.

for (m in c("observations", "metabarcoding")) {
  df_m <- df_all %>% filter(method == m)
  cat("\n==============================\n")
  cat("Method:", m, "\n")
  cat("==============================\n")

  # 1. Per-location network size
  cat("\n--- Network size per location ---\n")
  df_m %>%
    group_by(location) %>%
    summarise(
      n_plants      = n_distinct(lower_level),
      n_pollinators = n_distinct(higher_level),
      n_obs_links   = sum(ground_truth == 1),
      n_total_pairs = n(),
      .groups = "drop"
    ) %>% print()

  # 2. Species overlap across locations
  cat("\n--- Plant species shared across locations ---\n")
  plants_by_loc <- df_m %>%
    group_by(location) %>%
    summarise(plants = list(unique(lower_level)), .groups = "drop")
  locs <- plants_by_loc$location
  for (i in seq_along(locs)) {
    for (j in seq_along(locs)) {
      if (j > i) {
        shared <- length(intersect(plants_by_loc$plants[[i]],
                                   plants_by_loc$plants[[j]]))
        cat(sprintf("  %s ∩ %s: %d shared plants\n", locs[i], locs[j], shared))
      }
    }
  }

  cat("\n--- Pollinator species shared across locations ---\n")
  polls_by_loc <- df_m %>%
    group_by(location) %>%
    summarise(polls = list(unique(higher_level)), .groups = "drop")
  for (i in seq_along(locs)) {
    for (j in seq_along(locs)) {
      if (j > i) {
        shared <- length(intersect(polls_by_loc$polls[[i]],
                                   polls_by_loc$polls[[j]]))
        cat(sprintf("  %s ∩ %s: %d shared pollinators\n", locs[i], locs[j], shared))
      }
    }
  }

  # 3. Interaction ID overlap (all pairs, regardless of ground_truth)
  cat("\n--- Interaction ID overlap across locations (all pairs) ---\n")
  ids_by_loc <- df_m %>%
    group_by(location) %>%
    summarise(ids = list(unique(interaction_id)), .groups = "drop")
  for (i in seq_along(locs)) {
    for (j in seq_along(locs)) {
      if (j > i) {
        shared <- length(intersect(ids_by_loc$ids[[i]], ids_by_loc$ids[[j]]))
        cat(sprintf("  %s ∩ %s: %d shared pairs\n", locs[i], locs[j], shared))
      }
    }
  }

  # 4. n_obs_total distribution: how many interactions are seen in >1 location?
  cat("\n--- Distribution of n_obs_total (# locations where ground_truth == 1) ---\n")
  obs_counts <- df_m %>%
    group_by(interaction_id) %>%
    summarise(n_obs_total = sum(ground_truth == 1), .groups = "drop")
  print(as.data.frame(table(n_obs_total = obs_counts$n_obs_total)))

  # 5. Sample 5 observed interaction IDs from each location
  cat("\n--- Sample observed interaction IDs per location ---\n")
  df_m %>%
    filter(ground_truth == 1) %>%
    group_by(location) %>%
    slice_head(n = 5) %>%
    select(location, interaction_id) %>%
    print(n = Inf)
}

# ---- community composition plots ----
# Tile plot: species (y) × location (x), coloured by how many locations the species appears in.
# Plants and pollinators in separate facets; species ordered most-widespread → unique, then A–Z.

plot_community <- function(df_all, method_name) {
  df_m      <- df_all %>% filter(method == method_name)
  locations <- sort(unique(df_m$location))

  sp_present <- bind_rows(
    df_m %>% distinct(location, species = lower_level)  %>% mutate(taxon = "Plants"),
    df_m %>% distinct(location, species = higher_level) %>% mutate(taxon = "Pollinators")
  )

  n_locs_per_sp <- sp_present %>%
    group_by(taxon, species) %>%
    summarise(n_locs = n_distinct(location), .groups = "drop")

  # full grid with absences filled in
  plot_df <- sp_present %>%
    distinct(taxon, species) %>%
    tidyr::crossing(location = locations) %>%
    left_join(sp_present %>% mutate(present = TRUE),
              by = c("taxon", "species", "location")) %>%
    left_join(n_locs_per_sp, by = c("taxon", "species")) %>%
    mutate(
      present  = !is.na(present),
      fill_val = factor(
        ifelse(present, as.character(n_locs), "Absent"),
        levels = c("Absent", "1", "2", "3")
      )
    )

  # species ordered: most widespread first, then A–Z within each n_locs group
  sp_order <- n_locs_per_sp %>%
    arrange(taxon, desc(n_locs), species) %>%
    pull(species)

  plot_df <- plot_df %>%
    mutate(species = factor(species, levels = rev(sp_order)))

  n_species <- n_distinct(plot_df$species)

  ggplot(plot_df, aes(x = location, y = species, fill = fill_val)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_manual(
      values = c("Absent" = "grey92", "1" = "#aec6cf", "2" = "#f4a261", "3" = "#2a9d8f"),
      labels = c("Absent", "1 location", "2 locations", "3 locations"),
      name   = "Presence"
    ) +
    facet_wrap(~ taxon, scales = "free_y", ncol = 2) +
    labs(title = paste("Community composition —", method_name),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y      = element_text(size = max(4, min(8, 400 / n_species)), face = "italic"),
      axis.text.x      = element_text(size = 11, face = "bold"),
      strip.text       = element_text(size = 12, face = "bold"),
      panel.grid       = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_text(size = 9),
      plot.title       = element_text(face = "bold")
    )
}

gg_comm_obs <- plot_community(df_all, "observations")
gg_comm_obs

gg_comm_mb <- plot_community(df_all, "metabarcoding")
gg_comm_mb

# ---- interaction richness column plots ----
# One plot per method × location × taxon (6 plots per method).
# Bar length = number of observed interaction partners (ground_truth == 1).
# Species sorted by richness; all names shown on the y-axis.

make_richness_plot <- function(df_all, method_name, location_name, taxon_name) {
  obs <- df_all %>%
    filter(method == method_name, location == location_name, ground_truth == 1)

  richness <- if (taxon_name == "Plants") {
    obs %>%
      group_by(species = lower_level) %>%
      summarise(n_partners = n_distinct(higher_level), .groups = "drop")
  } else {
    obs %>%
      group_by(species = higher_level) %>%
      summarise(n_partners = n_distinct(lower_level), .groups = "drop")
  }

  richness <- richness %>%
    mutate(species = factor(species, levels = sort(unique(species), decreasing = TRUE)))

  ggplot(richness, aes(x = n_partners, y = species)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n_partners), hjust = -0.3, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title    = paste(location_name, "—", taxon_name),
      subtitle = method_name,
      x        = "Number of interaction partners",
      y        = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y        = element_text(size = 8, face = "italic"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(face = "bold"),
      plot.subtitle      = element_text(color = "grey40")
    )
}

# generate all plots and store in named lists
locations <- c("Starkey", "Threemile", "Zumwalt")
taxa      <- c("Plants", "Pollinators")
methods   <- c("observations", "metabarcoding")

richness_plots <- lapply(methods, function(m) {
  lapply(locations, function(loc) {
    lapply(taxa, function(tx) {
      make_richness_plot(df_all, m, loc, tx)
    }) %>% setNames(taxa)
  }) %>% setNames(locations)
}) %>% setNames(methods)

# display: one combined plot per method × taxon (locations side by side)
for (m in methods) {
  for (tx in taxa) {
    combined <- richness_plots[[m]][[locations[1]]][[tx]] |
                richness_plots[[m]][[locations[2]]][[tx]] |
                richness_plots[[m]][[locations[3]]][[tx]]
    print(combined + plot_annotation(title = paste(tx, "—", m),
                                     theme = theme(plot.title = element_text(face = "bold",
                                                                             size = 13))))
  }
}

# ---- classify and corroborate within one method ----
# For a given method, we seek for spatial evidence: an interaction is "Observed elsewhere"
# if ground_truth == 1 in at least one OTHER location under the SAME method.
# Each location's data is treated as an independent row throughout; counts are only
# summed for the proportional display in the alluvial plot (not for classification).

classify_by_location <- function(df_all, method_name) {

  df_method <- df_all %>% filter(method == method_name)   # fix: == not =

  # count how many locations (under this method) observed each interaction
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
      # obs_elsewhere: observed in ≥1 OTHER location regardless of this one
      # subtracting as.integer(original_binary == 1) removes this location's contribution
      obs_elsewhere = (n_obs_total - as.integer(original_binary == 1)) >= 1,

      validation = ifelse(obs_elsewhere, "Observed elsewhere", "Not observed elsewhere")
    )

  df_categorized <- df_flagged %>%
    mutate(
      link_category = case_when(
        # observed in only this location
        is_unique     & original_binary == 1 & predicted_bin == 1 ~ "locally_unique_links",
        is_unique     & original_binary == 1 & predicted_bin == 0 ~ "unsupported_links",

        # observed in ≥2 locations (including this one)
        is_shared     & original_binary == 1 & predicted_bin == 1 ~ "recurrent_links",
        is_shared     & original_binary == 1 & predicted_bin == 0 ~ "model_elusive_links",

        # never observed at any location
        is_all_zero   & original_binary == 0 & predicted_bin == 0 ~ "likely_forbidden",
        is_all_zero   & original_binary == 0 & predicted_bin == 1 ~ "unconfirmed_links",

        # absent here but observed elsewhere
        obs_elsewhere & original_binary == 0 & predicted_bin == 0 ~ "locally_absent_links",
        obs_elsewhere & original_binary == 0 & predicted_bin == 1 ~ "possibly_missing_links",

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

  # one row per (location × interaction × category) — no averaging
  summary_tbl <- df_categorized %>%
    filter(link_category != "unclassified", confusion != "UNK") %>%
    count(location, confusion, validation, link_category, name = "value")

  list(categorized = df_categorized, summary = summary_tbl)
}

# ---- run for both methods ----
results_obs <- classify_by_location(df_all, "observations")
results_mb  <- classify_by_location(df_all, "metabarcoding")

df_obs_categorized <- results_obs$categorized
summary_obs        <- results_obs$summary

df_mb_categorized  <- results_mb$categorized
summary_mb         <- results_mb$summary

# ---- alluvial plot aesthetics (shared) ----
col_confusion  <- c("TP" = "lightsteelblue",  "FP" = "lightsteelblue2",
                    "TN" = "rosybrown",        "FN" = "rosybrown2")
col_validation <- c("Observed elsewhere"     = "sandybrown",
                    "Not observed elsewhere" = "thistle3")
col_subcats    <- c(
  "recurrent_links"        = "coral3",
  "locally_unique_links"   = "thistle3",
  "model_elusive_links"          = "coral",
  "unsupported_links"      = "thistle1",
  "possibly_missing_links" = "coral2",
  "unconfirmed_links"         = "thistle",
  "locally_absent_links"         = "coral1",
  "likely_forbidden"       = "thistle2"
)
stratum_fill <- c(col_confusion, col_validation, col_subcats)

order_confusion  <- c("TP", "FP", "TN", "FN")
order_validation <- c("Not observed elsewhere", "Observed elsewhere")
order_subtypes   <- c(
  "locally_unique_links",   "unconfirmed_links",
  "likely_forbidden",       "unsupported_links",
  "recurrent_links",        "possibly_missing_links",
  "locally_absent_links",         "model_elusive_links"
)

metric          <- "prop"   # "value" for raw counts, "prop" for proportions
flow_alpha      <- 0.7
flow_colour     <- NA
bg_col          <- "white"
name_size       <- 4.2
label_nudge_x   <- 0.03
label_color_map <- c("Confusion"  = "mistyrose4",
                     "Validation" = "mistyrose4",
                     "Subtype"    = "mistyrose4")

# ---- alluvial plot function ----
# summary_tbl: output of classify_by_location()$summary (has location column)
# method_label: string used in plot title
# Counts are summed across locations for display only — classification was done per location.

make_alluvial <- function(summary_tbl, method_label,
                          order_val  = order_validation,
                          col_val    = col_validation,
                          val_subtitle = "Spatial corroboration") {

  local_fill <- c(col_confusion, col_val, col_subcats)

  # sum counts across locations: each (confusion × validation × category) cell
  # accumulates observations from all 3 locations; no averaging occurs
  flows_tbl <- summary_tbl %>%
    group_by(confusion, validation, link_category) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    mutate(
      L1          = confusion,
      L2          = validation,
      L3          = link_category,
      total       = sum(value),
      prop        = ifelse(total > 0, value / total, 0),
      alluvium_id = paste(L1, L3, sep = "⟂")
    )

  flows_long <- flows_tbl %>%
    select(L1, L2, L3, value, prop, alluvium_id) %>%
    tidyr::pivot_longer(c(L1, L2, L3), names_to = "axis", values_to = "stratum") %>%
    dplyr::mutate(
      axis    = dplyr::recode(axis, L1 = "Confusion", L2 = "Validation", L3 = "Subtype"),
      stratum = dplyr::case_when(
        axis == "Confusion"  ~ factor(stratum, levels = order_confusion),
        axis == "Validation" ~ factor(stratum, levels = order_val),
        axis == "Subtype"    ~ factor(stratum, levels = order_subtypes),
        TRUE ~ factor(stratum)
      ),
      axis = factor(axis, levels = c("Confusion", "Validation", "Subtype"))
    )

  yr <- diff(range(flows_long[[metric]], na.rm = TRUE))

  # stratum totals for labelling
  stratum_totals <- flows_long %>%
    dplyr::group_by(axis, stratum) %>%
    dplyr::summarise(val = sum(.data[[metric]]), .groups = "drop") %>%
    dplyr::mutate(stratum_chr = as.character(stratum))

  # hidden plot to extract stratum box geometry for label placement
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
    scale_fill_manual(values = local_fill, guide = "none") +
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
    scale_color_manual(values = label_color_map, guide = "none") +
    coord_cartesian(clip = "off") +
    labs(
      title    = paste("Link classification —", method_label),
      subtitle = paste("Confusion →", val_subtitle, "→ Link category"),
      caption  = "Proportions reflect per-location assessments — each interaction is assessed once per site.",
      x        = NULL,
      y        = if (metric == "prop") "Proportion of interactions" else "Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(size = 12, face = "bold"),
      plot.title         = element_text(face = "bold"),
      plot.subtitle      = element_text(color = "grey30"),
      plot.caption       = element_text(color = "grey50", size = 8, hjust = 0),
      plot.margin        = margin(20, 120, 20, 20)
    )
}

# ---- produce one alluvial plot per method ----
gg_obs <- make_alluvial(summary_obs, "Observations")
gg_obs

gg_mb <- make_alluvial(summary_mb, "Metabarcoding (LDB)")
gg_mb

# ============================================================================
# ---- FrugInt MN_2024 — spatial corroboration across two Pistacia sites ----
# Two mist-netting networks (Hato Ratón and Southern site) are used for
# spatial corroboration: an interaction is "Observed elsewhere" if it was
# recorded at the other Pistacia site under the same method (mist-netting).
# ============================================================================

# ---- FrugInt MN_2024 — load prediction data ----
frugint_HatoRaton <- read_csv(
  "results/predictions/frugint_MN2024_HatoRaton_prediction_results.csv",
  show_col_types = FALSE
) %>% mutate(location = "HatoRaton", method = "mist-netting")

frugint_South <- read_csv(
  "results/predictions/frugint_MN2024_South_prediction_results.csv",
  show_col_types = FALSE
) %>% mutate(location = "South", method = "mist-netting")

df_frugint <- bind_rows(frugint_HatoRaton, frugint_South) %>%
  mutate(interaction_id = paste(higher_level, lower_level, sep = "___"))

# ---- FrugInt MN_2024 — diagnostic checks ----
frugint_locs <- sort(unique(df_frugint$location))

cat("\n==============================\n")
cat("FrugInt MN_2024 — mist-netting\n")
cat("==============================\n")

cat("\n--- Network size per location ---\n")
df_frugint %>%
  group_by(location) %>%
  summarise(
    n_plants   = n_distinct(lower_level),
    n_animals  = n_distinct(higher_level),
    n_obs_links = sum(ground_truth == 1),
    n_total_pairs = n(),
    .groups = "drop"
  ) %>% print()

cat("\n--- Plant species shared across locations ---\n")
plants_fl <- df_frugint %>%
  group_by(location) %>% summarise(sp = list(unique(lower_level)), .groups = "drop")
for (i in seq_along(frugint_locs)) for (j in seq_along(frugint_locs)) if (j > i) {
  cat(sprintf("  %s ∩ %s: %d shared plants\n", frugint_locs[i], frugint_locs[j],
              length(intersect(plants_fl$sp[[i]], plants_fl$sp[[j]]))))
}

cat("\n--- Animal species shared across locations ---\n")
animals_fl <- df_frugint %>%
  group_by(location) %>% summarise(sp = list(unique(higher_level)), .groups = "drop")
for (i in seq_along(frugint_locs)) for (j in seq_along(frugint_locs)) if (j > i) {
  cat(sprintf("  %s ∩ %s: %d shared animals\n", frugint_locs[i], frugint_locs[j],
              length(intersect(animals_fl$sp[[i]], animals_fl$sp[[j]]))))
}

cat("\n--- Interaction ID overlap across locations (all pairs) ---\n")
ids_fl <- df_frugint %>%
  group_by(location) %>% summarise(ids = list(unique(interaction_id)), .groups = "drop")
for (i in seq_along(frugint_locs)) for (j in seq_along(frugint_locs)) if (j > i) {
  cat(sprintf("  %s ∩ %s: %d shared pairs\n", frugint_locs[i], frugint_locs[j],
              length(intersect(ids_fl$ids[[i]], ids_fl$ids[[j]]))))
}

cat("\n--- Distribution of n_obs_total (# locations where ground_truth == 1) ---\n")
df_frugint %>%
  group_by(interaction_id) %>%
  summarise(n_obs_total = sum(ground_truth == 1), .groups = "drop") %>%
  count(n_obs_total) %>% print()

# ---- FrugInt MN_2024 — community composition tile plot ----
gg_frugint_comm <- plot_community(df_frugint, "mist-netting")
gg_frugint_comm

# ---- FrugInt MN_2024 — interaction richness plots ----
frugint_taxa <- c("Plants", "Pollinators")
for (tx in frugint_taxa) {
  combined <- make_richness_plot(df_frugint, "mist-netting", frugint_locs[1], tx) |
              make_richness_plot(df_frugint, "mist-netting", frugint_locs[2], tx)
  print(combined + plot_annotation(
    title = paste("FrugInt MN_2024 —", tx),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  ))
}

# ---- FrugInt MN_2024 — classify and corroborate ----
results_frugint    <- classify_by_location(df_frugint, "mist-netting")
df_frugint_cat     <- results_frugint$categorized
summary_frugint    <- results_frugint$summary

# ---- FrugInt MN_2024 — alluvial plot ----
gg_frugint <- make_alluvial(summary_frugint, "MN_2024 mist-netting (Hato Ratón + South)")
gg_frugint

# ============================================================================
# ---- FrugInt — additional-method validation of spurious and forbidden links ----
# Extends the 3-axis alluvial with a 4th column (additional method evidence).
# For unconfirmed_links and likely_forbidden only, links are split into:
#   "confirmed"   — interaction observed (ground_truth == 1) in the additional method
#   "unconfirmed" — no evidence from the additional method
# All other subtypes flow to a single "Other categories" stratum on axis 4.
# The function is generic: pass any set of observed interaction_ids and a label.
# ============================================================================

# ---- FrugInt — load additional-method observed interactions (BC seed matrix) ----
bc_seed_mat <- read.csv(
  "data/adjacency_matrices/frugint_BCseed_Pistacia.csv",
  row.names = 1, check.names = FALSE
)

# rows = plants (lower_level), columns = animals (higher_level)
# interaction_id format matches df_frugint: paste(higher_level, lower_level, sep = "___")
add_method_obs_ids <- bc_seed_mat %>%
  tibble::rownames_to_column("plant") %>%
  tidyr::pivot_longer(-plant, names_to = "animal", values_to = "observed") %>%
  dplyr::filter(observed == 1) %>%
  dplyr::mutate(interaction_id = paste(animal, plant, sep = "___")) %>%
  dplyr::pull(interaction_id)

cat(sprintf("\nAdditional method (BC seed Pistacia matrix) — observed interactions: %d\n",
            length(add_method_obs_ids)))

# ---- FrugInt — 4-axis alluvial function with additional-method validation ----
# df_categorized  : output of classify_by_location()$categorized
# add_obs_ids     : character vector of interaction_ids observed in the additional method
# method_label    : string used in the plot title (primary method)
# validation_label: name of the additional method shown on axis 4

make_alluvial_validated <- function(df_categorized, add_obs_ids, method_label,
                                    validation_label = "Additional method") {

  # axis 4: one confirmed + one unconfirmed stratum per category, interleaved
  # so each axis-3 flow visually splits within its own category band.
  # Confirmed strata inherit the category's axis-3 color; unconfirmed = grey88.
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
      title    = paste("Link classification + additional method validation —",
                       method_label),
      subtitle = paste("Confusion → Validation → Subtype →", validation_label),
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

# ---- FrugInt — produce validated alluvial plot ----
gg_frugint_validated <- make_alluvial_validated(
  df_frugint_cat,
  add_method_obs_ids,
  method_label      = "MN_2024 mist-netting (Hato Ratón + South)",
  validation_label  = "BC seed (Pistacia)"
)
gg_frugint_validated

# ---- FrugInt MN_2024 — evidence summary across all categories ----
# has_evidence = spatial evidence (obs_elsewhere) OR orthogonal evidence (BC seed matrix)
df_frugint_cat %>%
  filter(link_category != "unclassified", confusion != "UNK") %>%
  mutate(has_evidence = obs_elsewhere | (interaction_id %in% add_method_obs_ids)) %>%
  count(link_category, has_evidence, name = "n") %>%
  pivot_wider(names_from = has_evidence, values_from = n,
              values_fill = 0, names_prefix = "ev_") %>%
  rename(with_evidence = ev_TRUE, without_evidence = ev_FALSE) %>%
  mutate(total    = with_evidence + without_evidence,
         pct_with = scales::percent(with_evidence / total, accuracy = 1)) %>%
  arrange(desc(with_evidence)) %>%
  print()

# ============================================================================
# ---- FrugInt metaweb — method corroboration with BC seed matrix ----
# The metaweb combines both MN_2024 Pistacia sites (Hato Ratón + South).
# Predictions are corroborated by checking whether each interaction was
# observed in the BC seed matrix (orthogonal sampling method).
# There is no second spatial replicate here; corroboration is method-only.
# ============================================================================

# ---- FrugInt metaweb — load prediction data ----
frugint_metaweb <- read_csv(
  "results/predictions/frugint_MN2024_metaweb_prediction_results.csv",
  show_col_types = FALSE
) %>%
  mutate(interaction_id = paste(higher_level, lower_level, sep = "___"))

# ---- FrugInt metaweb — classify by method corroboration ----
# Analogous to classify_by_location() but corroboration comes from an
# additional method rather than spatial replication.
classify_by_method <- function(df_single, add_obs_ids,
                                obs_label     = "Observed in add. method",
                                not_obs_label = "Not observed in add. method") {
  df_flagged <- df_single %>%
    mutate(
      original_binary = ground_truth,
      predicted_bin   = classification,
      obs_in_method   = interaction_id %in% add_obs_ids,
      validation      = ifelse(obs_in_method, obs_label, not_obs_label)
    )

  df_categorized <- df_flagged %>%
    mutate(
      link_category = case_when(
        original_binary == 1 & predicted_bin == 1 &  obs_in_method ~ "recurrent_links",
        original_binary == 1 & predicted_bin == 1 & !obs_in_method ~ "locally_unique_links",
        original_binary == 1 & predicted_bin == 0 &  obs_in_method ~ "model_elusive_links",
        original_binary == 1 & predicted_bin == 0 & !obs_in_method ~ "unsupported_links",
        original_binary == 0 & predicted_bin == 1 &  obs_in_method ~ "possibly_missing_links",
        original_binary == 0 & predicted_bin == 1 & !obs_in_method ~ "unconfirmed_links",
        original_binary == 0 & predicted_bin == 0 &  obs_in_method ~ "locally_absent_links",
        original_binary == 0 & predicted_bin == 0 & !obs_in_method ~ "likely_forbidden",
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
    count(confusion, validation, link_category, name = "value")

  list(categorized = df_categorized, summary = summary_tbl)
}

results_metaweb <- classify_by_method(
  frugint_metaweb, add_method_obs_ids,
  obs_label     = "Observed in BC seed",
  not_obs_label = "Not observed in BC seed"
)
df_metaweb_cat  <- results_metaweb$categorized
summary_metaweb <- results_metaweb$summary

# ---- FrugInt metaweb — alluvial plot ----
order_val_bc <- c("Not observed in BC seed", "Observed in BC seed")
col_val_bc   <- c("Observed in BC seed"     = "sandybrown",
                  "Not observed in BC seed"  = "thistle3")

gg_metaweb_bc <- make_alluvial(
  summary_metaweb,
  method_label  = "MN_2024 metaweb (Hato Ratón + South)",
  order_val     = order_val_bc,
  col_val       = col_val_bc,
  val_subtitle  = "Method corroboration (BC seed)"
)
gg_metaweb_bc

# ---- FrugInt metaweb — BC seed evidence summary across all categories ----
# In the metaweb (single network), evidence comes only from the orthogonal method.
df_metaweb_cat %>%
  filter(link_category != "unclassified", confusion != "UNK") %>%
  mutate(has_evidence = obs_in_method) %>%
  count(link_category, has_evidence, name = "n") %>%
  pivot_wider(names_from = has_evidence, values_from = n,
              values_fill = 0, names_prefix = "ev_") %>%
  rename(with_evidence = ev_TRUE, without_evidence = ev_FALSE) %>%
  mutate(total    = with_evidence + without_evidence,
         pct_with = scales::percent(with_evidence / total, accuracy = 1)) %>%
  arrange(desc(with_evidence)) %>%
  print()

# ============================================================================
# ---- FrugInt — link category comparison: spatial vs metaweb analysis ----
# Spatial analysis: MN_2024 two-site spatial corroboration (counts are per
#   location-level rows, consistent with the alluvial plots above).
# Metaweb analysis: single combined network corroborated by BC seed matrix.
# Proportions are used for comparison because the two analyses have different
# denominators (2 × n_interactions vs 1 × n_interactions).
# ============================================================================

spatial_cats <- df_frugint_cat %>%
  filter(link_category != "unclassified", confusion != "UNK") %>%
  distinct(interaction_id, link_category) %>%   # unique interactions — df_frugint_cat has one
  count(link_category, name = "n") %>%          # row per location, so pairs at both sites
  mutate(analysis = "Spatial (2 sites)",        # would be double-counted without this
         prop     = n / sum(n))

metaweb_cats <- df_metaweb_cat %>%
  filter(link_category != "unclassified", confusion != "UNK") %>%
  count(link_category, name = "n") %>%
  mutate(analysis = "Metaweb (BC seed)",
         prop     = n / sum(n))

comparison_cats <- bind_rows(spatial_cats, metaweb_cats) %>%
  mutate(
    link_category  = factor(link_category, levels = order_subtypes),
    category_label = link_category %>%
      as.character() %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_to_sentence(),
    category_label = factor(category_label,
                            levels = stringr::str_replace_all(order_subtypes, "_", " ") %>%
                              stringr::str_to_sentence())
  )

# ---- comprehensive comparison table: unique interactions, all categories ----
# Spatial: distinct(interaction_id, confusion, validation, link_category) removes
#   per-location double-counting (recurrent/cryptic links appear at both sites).
# Metaweb: one row per interaction by construction.
# "corroboration" harmonises the validation axis: spatial "Observed elsewhere" and
#   metaweb "Observed in BC seed" both map to "With corroboration".
# bc_confirmed: for spurious/likely_forbidden only — whether the interaction was
#   observed in the BC seed matrix. In the metaweb analysis this is always
#   "unconfirmed" by construction (confirmed FPs become possibly_missing_links).

spatial_breakdown <- df_frugint_cat %>%
  filter(link_category != "unclassified", confusion != "UNK") %>%
  distinct(interaction_id, confusion, validation, link_category) %>%
  mutate(
    corroboration = ifelse(stringr::str_starts(validation, "Observed"),
                           "With corroboration", "Without corroboration"),
    bc_confirmed  = case_when(
      link_category %in% c("unconfirmed_links", "likely_forbidden") &
          (interaction_id %in% add_method_obs_ids)  ~ "confirmed",
      link_category %in% c("unconfirmed_links", "likely_forbidden") &
        !(interaction_id %in% add_method_obs_ids)   ~ "unconfirmed",
      TRUE ~ NA_character_
    ),
    analysis = "Spatial (2 sites)"
  )

metaweb_breakdown <- df_metaweb_cat %>%
  filter(link_category != "unclassified", confusion != "UNK") %>%
  mutate(
    corroboration = ifelse(obs_in_method, "With corroboration", "Without corroboration"),
    bc_confirmed  = case_when(
      link_category %in% c("unconfirmed_links", "likely_forbidden") &  obs_in_method  ~ "confirmed",
      link_category %in% c("unconfirmed_links", "likely_forbidden") & !obs_in_method  ~ "unconfirmed",
      TRUE ~ NA_character_
    ),
    analysis = "Metaweb (BC seed)"
  )

# Verification: rows vs unique interaction_ids for spatial spurious/forbidden
cat("\n--- Verification: spatial per-location rows vs unique interaction_ids ---\n")
df_frugint_cat %>%
  filter(link_category %in% c("unconfirmed_links", "likely_forbidden"), confusion != "UNK") %>%
  group_by(link_category) %>%
  summarise(
    n_rows          = n(),
    n_unique        = n_distinct(interaction_id),
    n_at_both_sites = n_rows - n_unique,
    .groups = "drop"
  ) %>% print()

# Verification: metaweb spurious/forbidden should have 0 BC-seed-confirmed (by construction)
cat("\n--- Verification: metaweb spurious/forbidden — BC seed confirmed count (expect 0) ---\n")
df_metaweb_cat %>%
  filter(link_category %in% c("unconfirmed_links", "likely_forbidden"), confusion != "UNK") %>%
  group_by(link_category) %>%
  summarise(n_confirmed = sum(obs_in_method), .groups = "drop") %>% print()

cat("\n--- Comprehensive comparison: unique interactions, all categories ---\n")
bind_rows(
  spatial_breakdown %>%
    select(interaction_id, analysis, confusion, corroboration, link_category, bc_confirmed),
  metaweb_breakdown %>%
    select(interaction_id, analysis, confusion, corroboration, link_category, bc_confirmed)
) %>%
  count(analysis, confusion, corroboration, link_category, bc_confirmed, name = "n") %>%
  group_by(analysis) %>%
  mutate(pct = scales::percent(n / sum(n), accuracy = 0.1)) %>%
  ungroup() %>%
  mutate(cell = paste0(n, " (", pct, ")")) %>%
  select(-n, -pct) %>%
  tidyr::pivot_wider(names_from = analysis, values_from = cell, values_fill = "0 (0%)") %>%
  arrange(
    confusion,
    factor(corroboration, levels = c("With corroboration", "Without corroboration")),
    link_category,
    bc_confirmed
  ) %>%
  print(n = Inf)

# ---- bar chart ----
ggplot(comparison_cats,
       aes(x = category_label, y = prop, fill = analysis)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = scales::percent(prop, accuracy = 1)),
            position = position_dodge(width = 0.75),
            vjust = -0.4, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = c("Spatial (2 sites)"  = "steelblue3",
                                "Metaweb (BC seed)" = "coral3")) +
  labs(
    title    = "Link category distribution: spatial vs metaweb analysis",
    subtitle = "Spatial: MN_2024 Hato Ratón + South | Metaweb: single combined network",
    x        = NULL,
    y        = "Proportion of interactions",
    fill     = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x        = element_text(angle = 30, hjust = 1, size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "top",
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(color = "grey30")
  )
