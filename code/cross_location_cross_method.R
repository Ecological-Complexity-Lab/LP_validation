# ---- cross location cross method analysis ----
# here we corroborate predictions produced by Barry's ecoILP application using data from the
# plant-pollinator networks in 3 locations that were samples using two different
# methods (observations, DNA metabarcoding)
# network data from Arstingstall, Katherine A., Sandra J. DeBano, Xiaoping Li, David E. Wooster, Mary M. Rowland, Skyler Burrows, and Kenneth Frost. 2021. "Capabilities and Limitations of Using DNA Metabarcoding to Study Plant-Pollinator Interactions." Molecular Ecology 30 (20): 5266–97.
# predictions done on the full matrix (no link withholding)
# classification threshold is 0.5

# ---- load libraries ----
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggalluvial)
library(scales)

# ---- read prediction data ----
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

# ---- classify and corroborate within one method ----
# For a given method, corroboration is purely spatial: an interaction is "Corroborated"
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

      validation = ifelse(obs_elsewhere, "Corroborated", "Not corroborated")
    )

  df_categorized <- df_flagged %>%
    mutate(
      link_category = case_when(
        # observed in only this location
        is_unique     & original_binary == 1 & predicted_bin == 1 ~ "locally_unique_links",
        is_unique     & original_binary == 1 & predicted_bin == 0 ~ "unsupported_links",

        # observed in ≥2 locations (including this one)
        is_shared     & original_binary == 1 & predicted_bin == 1 ~ "recurrent_links",
        is_shared     & original_binary == 1 & predicted_bin == 0 ~ "cryptic_links",

        # never observed at any location
        is_all_zero   & original_binary == 0 & predicted_bin == 0 ~ "likely_forbidden",
        is_all_zero   & original_binary == 0 & predicted_bin == 1 ~ "spurious_links",

        # absent here but observed elsewhere
        obs_elsewhere & original_binary == 0 & predicted_bin == 0 ~ "feasible_links",
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
col_validation <- c("Corroborated"     = "sandybrown",
                    "Not corroborated" = "thistle3")
col_subcats    <- c(
  "recurrent_links"        = "coral3",
  "locally_unique_links"   = "thistle3",
  "cryptic_links"          = "coral",
  "unsupported_links"      = "thistle1",
  "possibly_missing_links" = "coral2",
  "spurious_links"         = "thistle",
  "feasible_links"         = "coral1",
  "likely_forbidden"       = "thistle2"
)
stratum_fill <- c(col_confusion, col_validation, col_subcats)

order_confusion  <- c("TP", "FP", "TN", "FN")
order_validation <- c("Not corroborated", "Corroborated")
order_subtypes   <- c(
  "locally_unique_links",   "spurious_links",
  "likely_forbidden",       "unsupported_links",
  "recurrent_links",        "possibly_missing_links",
  "feasible_links",         "cryptic_links"
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

make_alluvial <- function(summary_tbl, method_label) {

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
        axis == "Validation" ~ factor(stratum, levels = order_validation),
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
    scale_fill_manual(values = stratum_fill, guide = "none") +
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
      subtitle = "Confusion → Spatial corroboration → Link category",
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
      plot.margin        = margin(20, 120, 20, 20)
    )
}

# ---- produce one alluvial plot per method ----
gg_obs <- make_alluvial(summary_obs, "Observations")
gg_obs

gg_mb <- make_alluvial(summary_mb, "Metabarcoding (LDB)")
gg_mb
