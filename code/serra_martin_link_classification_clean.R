# ---- Serra-Martin link classification and alluvial plots (reorganised) ----
# Loads cross-site SVD self-prediction results (6 habitat patches x 2 sampling
# methods), classifies every predicted link via spatial corroboration across
# the 6 sites, then corroborates each method's classifications using evidence
# from the other sampling method (obs vs. rpi camera).
#
# Outputs:
#   - 4 static Sankey/alluvial figures (labeled + unlabeled, one per method)
#   - 2 "all categories split" inspection Sankeys (split_all = TRUE)
#   - 1 interactive Plotly Sankey (threshold + method sliders)
#   - Tile maps for the most species-rich site (overall, and for the richest
#     camera/rpi site specifically)
#   - A text summary of the unobserved-link breakdown per method
#
# ---- 1. Libraries ----
# Checks that every package used below is installed and installs anything
# missing before loading it, so this script runs standalone on a fresh
# machine. ggsankey isn't on CRAN, so it's installed from GitHub via remotes.

required_cran_packages <- c(
  "dplyr", "readr", "tidyr", "stringr", "ggplot2",
  "plotly", "htmlwidgets", "jsonlite",         # used by interactive_sankey/make_sankey.R
  "remotes"                                    # needed to install ggsankey from GitHub
)

missing_cran <- setdiff(required_cran_packages, rownames(installed.packages()))
if (length(missing_cran) > 0) {
  message("Installing missing packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
}

if (!requireNamespace("ggsankey", quietly = TRUE)) {
  message("Installing ggsankey from GitHub (davidsjoberg/ggsankey)...")
  remotes::install_github("davidsjoberg/ggsankey")
}

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggsankey)

# ---- 2. Aesthetics ----

# -- Sankey node/flow colors, shared by classification (Confusion/Validation/
#    Subtype axes) and the maps below (which reuse the same subtype palette
#    via col_subcats through map_colors) --
col_confusion <- c("TP" = "lightsteelblue1",  "FP" = "lightsteelblue",
                   "TN" = "rosybrown",        "FN" = "rosybrown2")
col_validation <- c("Observed elsewhere"     = "sandybrown",
                    "Not observed elsewhere" = "thistle3")
col_subcats <- c(
  "recurrent"        = "coral",
  "locally_unique"   = "thistle3",
  "model_elusive"    = "coral2",
  "weakly-supported"      = "thistle",
  "possibly_missing" = "coral1",
  "phantom"          = "thistle",
  "locally_absent"   = "coral3",
  "possibly_forbidden" = "thistle2"
)

flow_alpha <- 0.7  # Sankey flow (ribbon) transparency

# -- Tile-map colors: the original map palette, regrouped onto the contextual
#    evidence axis rather than the confusion quadrant.
#      warm (corals, rose, amber) = observed in >=1 replicate
#      cool (teals, blue, purple) = not observed in any replicate
#      aquamarine                 = corroborated by the other sampling method
#    Within warm, dark marks links also recorded at this site; within cool, the
#    teals mark links recorded here and the blue/purple those recorded nowhere.
#    The two aquamarines take the place of the palette's former light-coral and
#    light-lavender have-evidence tints; the other eight colors are unchanged,
#    only reassigned. Warm/cool membership was checked against the validation axis. --
map_display_levels <- c(
  "Locally unique",
  "Recurrent",
  "Phantom — have evidence", "Phantom — no evidence",
  "Possibly missing",
  "Weakly-supported",
  "Model elusive",
  "Possibly forbidden — have evidence", "Possibly forbidden — no evidence",
  "Locally absent"
)

map_colors <- c(
  # -- Observed elsewhere (recorded in >=1 replicate) — warm --
  #    dark = also recorded here, light = recorded only elsewhere
  "Recurrent"                          = "#C05030",  # dark coral   also recorded here, predicted
  "Model elusive"                      = "#8A3050",  # deep rose    also recorded here, not predicted
  "Possibly missing"                   = "#DBA040",  # amber        only elsewhere, predicted
  "Locally absent"                     = "#F9C8B0",  # lightsteelblue   only elsewhere, not predicted
  # -- Not observed elsewhere — cool --
  #    teal = recorded here, blue/purple = recorded nowhere
  "Locally unique"                     = "#B0C4DE",  # light teal   recorded here only, predicted
  "Weakly-supported"                   = "#E8A8B8",  #    recorded here only, not predicted
  "Phantom — no evidence"              = "#D5C8F0",  # recorded nowhere, predicted
  "Possibly forbidden — no evidence"   = "#5E3DA0",  # dark purple  recorded nowhere, not predicted
  # -- Corroborated by the other sampling method — aquamarine (as in the Sankey) --
  #    these replace the palette's two former have-evidence tints
  "Phantom — have evidence"            = "#76C7C0",   # #7FFFD4
  "Possibly forbidden — have evidence" = "#1F7E72"   # #66CDAA
)


# -- Lavender gradient for the unobserved subtypes (locally_unique lightest ->
#    possibly_forbidden darkest). Defined here rather than inside
#    make_sankey_validated() so the Sankey and the Sankey-colored map below read
#    from one definition and cannot drift apart. --
col_lav <- c("locally_unique"     = "#C3BAD5",
             "phantom"            = "#AFA2C4",
             "weakly-supported"   = "#9B8BB4",
             "possibly_forbidden" = "#8878A4")

# -- Alternative tile-map palette: the alluvial (Sankey) colors, so each map cell
#    carries the same color as its ribbon in sankey_obs_validated_with_cameras.
#    Same 10 display levels as map_colors, and assembled from the very objects
#    make_sankey_validated() uses for axis 4 (col_lav, col_subcats, and the
#    aquamarine have-evidence pair). Used by map_rich_sankey in Section 9. --
map_colors_sankey <- c(
  "Locally unique"                     = unname(col_lav["locally_unique"]),
  "Recurrent"                          = unname(col_subcats["recurrent"]),
  "Phantom — have evidence"            = "aquamarine",
  "Phantom — no evidence"              = unname(col_lav["phantom"]),
  "Possibly missing"                   = unname(col_subcats["possibly_missing"]),
  "Weakly-supported"                   = unname(col_lav["weakly-supported"]),
  "Model elusive"                      = unname(col_subcats["model_elusive"]),
  "Possibly forbidden — have evidence"  = "aquamarine3",
  "Possibly forbidden — no evidence"    = unname(col_lav["possibly_forbidden"]),
  "Locally absent"                     = unname(col_subcats["locally_absent"])
)

# ---- 3. Functions ----

# Classifies every link prediction for one method (obs or rpi) using spatial
# corroboration: how many of the 6 sites recorded this species pair (ground
# truth), regardless of what happens at the focal site.
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
        is_unique     & original_binary == 1 & predicted_bin == 0 ~ "weakly-supported",
        is_shared     & original_binary == 1 & predicted_bin == 1 ~ "recurrent",
        is_shared     & original_binary == 1 & predicted_bin == 0 ~ "model_elusive",
        is_all_zero   & original_binary == 0 & predicted_bin == 0 ~ "possibly_forbidden",
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

# Builds one 4-axis Sankey (Within-network evaluation -> Contextual evidence ->
# Link category -> validation by the other method) for a single sampling
# method. By default (split_all = FALSE), only "phantom", "weakly-supported" and
# "possibly_forbidden" are split into have-evidence / no-evidence at axis 4 —
# those are the only subtypes ambiguous enough to benefit from cross-method
# corroboration. split_all = TRUE splits all 8 subtypes instead, for
# inspecting how corroboration lines up with every category at once.
make_sankey_validated <- function(df_categorized, add_obs_ids, method_label,
                                   validation_label = "Additional method",
                                   show_labels = TRUE,
                                   split_all = FALSE) {

  # ---- lavender gradient (locally_unique lightest -> possibly_forbidden darkest) --
  # Defined in Section 2 so map_colors_sankey reuses exactly these colors.
  lav <- col_lav

  # ---- axis 4 ordering, colors, node IDs, and rank vectors -------------------
  # split_all = FALSE (default): only phantom, weakly-supported, possibly_forbidden
  #   are split into have / no evidence; all other categories pass through.
  # split_all = TRUE  (inspection): all 8 categories are split.
  if (!split_all) {

    order_L4 <- c(
      "locally_unique",
      "recurrent",
      "phantom — have evidence",         "phantom — no evidence",
      "possibly_missing",
      "weakly-supported",
      "model_elusive",
      "possibly_forbidden — have evidence","possibly_forbidden — no evidence",
      "locally_absent"
    )

    col_L4 <- c(
      "locally_unique"                    = unname(lav["locally_unique"]),
      "recurrent"                         = unname(col_subcats["recurrent"]),
      "phantom — have evidence"           = "aquamarine",
      "phantom — no evidence"             = unname(lav["phantom"]),
      "possibly_missing"                  = unname(col_subcats["possibly_missing"]),
      "weakly-supported"                      = unname(lav["weakly-supported"]),
      "model_elusive"                     = unname(col_subcats["model_elusive"]),
      "possibly_forbidden — have evidence"  = "aquamarine3",
      "possibly_forbidden — no evidence"    = unname(lav["possibly_forbidden"]),
      "locally_absent"                    = unname(col_subcats["locally_absent"])
    )

    node_id_vec <- c(
      TN = 1L, FN = 2L, FP = 3L, TP = 4L,
      "Not observed elsewhere" = 5L, "Observed elsewhere" = 6L,
      "possibly_forbidden"                  =  7L,
      "possibly_forbidden — no evidence"    =  8L,
      "possibly_forbidden — have evidence"  =  9L,
      "weakly-supported"                      = 10L,
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
      possibly_forbidden = 10,
        "possibly_forbidden — no evidence" = 11, "possibly_forbidden — have evidence" = 12,
      "weakly-supported" = 20,
      phantom = 30,
        "phantom — no evidence" = 31, "phantom — have evidence" = 32,
      locally_unique = 40, locally_absent = 50,
      model_elusive = 60, possibly_missing = 70, recurrent = 80
    )

    node_colors_l4 <- c(
      "phantom — have evidence"           = unname(col_L4["phantom — have evidence"]),
      "phantom — no evidence"             = unname(col_L4["phantom — no evidence"]),
      "possibly_forbidden — have evidence"  = unname(col_L4["possibly_forbidden — have evidence"]),
      "possibly_forbidden — no evidence"    = unname(col_L4["possibly_forbidden — no evidence"])
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
      "weakly-supported — have evidence",      "weakly-supported — no evidence",
      "model_elusive — have evidence",     "model_elusive — no evidence",
      "possibly_forbidden — have evidence",  "possibly_forbidden — no evidence",
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
      "weakly-supported — have evidence"      = hev,
      "weakly-supported — no evidence"        = unname(lav["weakly-supported"]),
      "model_elusive — have evidence"     = hev,
      "model_elusive — no evidence"       = unname(col_subcats["model_elusive"]),
      "possibly_forbidden — have evidence"  = hev,
      "possibly_forbidden — no evidence"    = unname(lav["possibly_forbidden"]),
      "locally_absent — have evidence"    = hev,
      "locally_absent — no evidence"      = unname(col_subcats["locally_absent"])
    )

    # Axis-3 node IDs kept as-is (7–20); axis-4 split IDs start at 21.
    node_id_vec <- c(
      TN = 1L, FN = 2L, FP = 3L, TP = 4L,
      "Not observed elsewhere" = 5L, "Observed elsewhere" = 6L,
      "possibly_forbidden"                  =  7L,
      "weakly-supported"                      = 10L,
      "phantom"                           = 13L,
      "locally_unique"                    = 16L,
      "locally_absent"                    = 17L,
      "model_elusive"                     = 18L,
      "possibly_missing"                  = 19L,
      "recurrent"                         = 20L,
      "possibly_forbidden — no evidence"    = 21L,
      "possibly_forbidden — have evidence"  = 22L,
      "weakly-supported — no evidence"        = 23L,
      "weakly-supported — have evidence"      = 24L,
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
      possibly_forbidden = 10,
        "possibly_forbidden — no evidence" = 11, "possibly_forbidden — have evidence" = 12,
      "weakly-supported" = 20,
        "weakly-supported — no evidence" = 21, "weakly-supported — have evidence" = 22,
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
      "weakly-supported — have evidence"      = hev,
      "weakly-supported — no evidence"        = unname(lav["weakly-supported"]),
      "model_elusive — have evidence"     = hev,
      "model_elusive — no evidence"       = unname(col_subcats["model_elusive"]),
      "possibly_forbidden — have evidence"  = hev,
      "possibly_forbidden — no evidence"    = unname(lav["possibly_forbidden"]),
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
          link_category == "possibly_forbidden" & add_obs  ~ "possibly_forbidden — have evidence",
          link_category == "possibly_forbidden" & !add_obs ~ "possibly_forbidden — no evidence",
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
          link_category == "weakly-supported"     & add_obs  ~ "weakly-supported — have evidence",
          link_category == "weakly-supported"     & !add_obs ~ "weakly-supported — no evidence",
          link_category == "model_elusive"    & add_obs  ~ "model_elusive — have evidence",
          link_category == "model_elusive"    & !add_obs ~ "model_elusive — no evidence",
          link_category == "possibly_forbidden" & add_obs  ~ "possibly_forbidden — have evidence",
          link_category == "possibly_forbidden" & !add_obs ~ "possibly_forbidden — no evidence",
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
    "weakly-supported"           = unname(lav["weakly-supported"]),
    "model_elusive"          = unname(col_subcats["model_elusive"]),
    "possibly_forbidden"       = unname(lav["possibly_forbidden"]),
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
    # only the have/no-evidence splits (marked with "—") get one at axis 4.
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
    theme_sankey(base_size = 12) +
    theme(
      axis.text.x   = if (show_labels) element_text(size = 11, face = "bold") else element_blank(),
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30"),
      plot.margin   = margin(20, 20, 20, 20)
    )
}

# Summarises the unobserved-link breakdown for one method, pooling all 6
# sites. See the section comment above the calls (Summary section) for the
# full terminology (link-site records vs. unique species interactions, and
# the FP/TN subtype definitions).
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
    "    possibly_forbidden = predicted absent, never recorded at any site under %s\n",
    method_name))
  cat(sprintf(
    "    locally_absent   = predicted absent, recorded at ≥1 other site under %s\n\n",
    method_name))
  print_tbl(tn_tbl, n_tn, "TN")
  cat("\n")

  invisible(list(fp = fp_tbl, tn = tn_tbl))
}

# ---- 4. Data ----
# Cross-site SVD self-prediction results: every interaction, per site (6
# habitat patches) and per method (obs / rpi), with the ground-truth label,
# predicted probability and binary classification at the chosen threshold.
df_all <- read_csv(
  "results/predictions/serra_martin_loo_prediction_results.csv",
  show_col_types = FALSE
) %>%
  rename(location = focal_site, classification = prediction) %>%
  mutate(interaction_id = paste(higher_level, lower_level, sep = "___"))

# ---- 5. Classifications ----

# Per-location classification (contextual evidence from the 6 sites),
# separately for each sampling method.
results_obs <- classify_by_location(df_all, "obs")
results_rpi <- classify_by_location(df_all, "rpi")

df_obs_categorized <- results_obs$categorized
summary_obs        <- results_obs$summary

df_rpi_categorized <- results_rpi$categorized
summary_rpi        <- results_rpi$summary

# Cross-method evidence: interaction IDs observed (ground_truth == 1) by the
# *other* method, used to flag "have evidence" for ambiguous subtypes below.
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

# ---- 6. Diagnostics (sanity checks) ----
# Confirms record counts, site/method coverage, and that no interaction was
# recorded at more than the maximum possible 6 sites, before any plotting.

sites   <- sort(unique(df_all$location))
methods <- sort(unique(df_all$method))

cat("\n=== DIAGNOSTICS: Serra-Martin link classification ===\n\n")
cat(sprintf("Sites (n = %d): %s\n", length(sites), paste(sites, collapse = ", ")))
cat(sprintf("Methods (n = %d): %s\n\n", length(methods), paste(methods, collapse = ", ")))

cat("--- Link-site records per site x method ---\n")
df_all %>%
  count(location, method, name = "n_records") %>%
  print(n = Inf)

cat("\n--- Totals per method (all sites pooled) ---\n")
df_all %>%
  group_by(method) %>%
  summarise(
    n_records      = n(),
    n_sites        = n_distinct(location),
    n_observed     = sum(ground_truth == 1),
    n_unique_pairs = n_distinct(interaction_id),
    .groups = "drop"
  ) %>%
  print()

cat("\n--- Distribution of n_obs_total (# sites where ground_truth == 1), per method ---\n")
for (m in methods) {
  cat(sprintf("\nMethod: %s\n", m))
  obs_counts_m <- df_all %>%
    filter(method == m) %>%
    group_by(interaction_id) %>%
    summarise(n_obs_total = sum(ground_truth == 1), .groups = "drop")
  print(as.data.frame(table(n_obs_total = obs_counts_m$n_obs_total)))
  max_sites_seen <- max(obs_counts_m$n_obs_total)
  cat(sprintf("  Max sites any interaction was observed at: %d (must be <= %d)\n",
              max_sites_seen, length(sites)))
  stopifnot(max_sites_seen <= length(sites))
}

cat("\n--- Classification integrity check (every row must be classified) ---\n")
for (nm in c("obs", "rpi")) {
  df_cat <- if (nm == "obs") df_obs_categorized else df_rpi_categorized
  n_unclassified <- sum(df_cat$link_category == "unclassified")
  n_unk          <- sum(df_cat$confusion == "UNK")
  cat(sprintf("  %s: %d rows, %d unclassified, %d UNK confusion (both should be 0)\n",
              nm, nrow(df_cat), n_unclassified, n_unk))
  stopifnot(n_unclassified == 0, n_unk == 0)
}

cat(sprintf(
  "\nDiagnostics passed: record counts check out, no interaction exceeds %d sites, and every row was classified.\n\n",
  length(sites)
))

# ---- 7. Alluvial plots (Sankey diagrams) ----
# The 4 figures actually used downstream: labeled + unlabeled ("clean")
# Sankeys for each sampling method, each validated against the other method.

gg_obs_sankey <- make_sankey_validated(
  df_obs_categorized,
  add_obs_ids      = add_obs_for_obs,
  method_label     = "Direct observation (obs)",
  validation_label = "Raspberry Pi camera (rpi)"
)
gg_obs_sankey

# ggsave("results/figures/sankey_obs_validated_with_cameras.pdf",
#        plot = gg_obs_sankey, width = 14, height = 6, device = cairo_pdf)
# ggsave("results/figures/sankey_obs_validated_with_cameras.png",
#        plot = gg_obs_sankey, width = 14, height = 6, dpi = 300, bg = "white")

gg_obs_sankey_clean <- make_sankey_validated(
  df_obs_categorized,
  add_obs_ids      = add_obs_for_obs,
  method_label     = "Direct observation (obs)",
  validation_label = "Raspberry Pi camera (rpi)",
  show_labels      = FALSE             # <-- toggle here
)
gg_obs_sankey_clean
# ggsave("results/figures/sankey_obs_validated_with_cameras_clean.pdf",
#        plot = gg_obs_sankey_clean, width = 14, height = 6, device = cairo_pdf)
# ggsave("results/figures/sankey_obs_validated_with_cameras_clean.png",
#        plot = gg_obs_sankey_clean, width = 14, height = 6, dpi = 300, bg = "white")

ggsave("results/figures/gg_obs_sankey_clean.svg",
       plot = gg_obs_sankey_clean, width = 14, height = 10, bg = "transparent")

gg_rpi_sankey <- make_sankey_validated(
  df_rpi_categorized,
  add_obs_ids      = add_obs_for_rpi,
  method_label     = "Raspberry Pi camera (rpi)",
  validation_label = "Direct observation (obs)"
)
gg_rpi_sankey
# ggsave("results/figures/sankey_rpi_validated_with_obs.pdf",
#        plot = gg_rpi_sankey, width = 14, height = 6, device = cairo_pdf)
# ggsave("results/figures/sankey_rpi_validated_with_obs.png",
#        plot = gg_rpi_sankey, width = 14, height = 6, dpi = 300, bg = "white")

gg_rpi_sankey_clean <- make_sankey_validated(
  df_rpi_categorized,
  add_obs_ids      = add_obs_for_rpi,
  method_label     = "Raspberry Pi camera (rpi)",
  validation_label = "Direct observation (obs)",
  show_labels      = FALSE             # <-- toggle here
)
gg_rpi_sankey_clean
# ggsave("results/figures/sankey_rpi_validated_with_obs_clean.pdf",
#        plot = gg_rpi_sankey_clean, width = 14, height = 6, device = cairo_pdf)
# ggsave("results/figures/sankey_rpi_validated_with_obs_clean.png",
#        plot = gg_rpi_sankey_clean, width = 14, height = 6, dpi = 300, bg = "white")

# -- Inspection versions: all 8 categories split into have / no evidence --
# split_all = TRUE adds have evidence / no evidence at axis 4 for every
# category, not just phantom / weakly-supported / possibly_forbidden. Used to
# check corroboration patterns across every subtype at once.
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

# ---- 8. Interactive alluvial plot ----
# Self-contained interactive Plotly Sankey with method/threshold controls.
# Sourced as a separate script (from the project root) because it precomputes
# its own link datasets across 11 thresholds x 2 methods.
# Requires: plotly, htmlwidgets, jsonlite (checked/installed in Section 1).
source("code/make_interactive_sankey_app.R")

# ---- 9. Maps ----

# -- Most species-rich site (used in the paper) --
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

# Species ordered by overall degree across all sites under this method, so
# generalists cluster together regardless of what happens at this one site.
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
      link_category == "weakly-supported"      ~ "Weakly-supported",
      link_category == "possibly_forbidden" & add_obs  ~ "Possibly forbidden — have evidence",
      link_category == "possibly_forbidden" & !add_obs ~ "Possibly forbidden — no evidence",
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

# -- Same map, alluvial (Sankey) colors --
# Identical data, species ordering and theme as map_rich; only the fill palette
# changes, so every cell carries the color its ribbon has in the Sankey. Built by
# swapping the scale on map_rich rather than rebuilding the plot, so the two stay
# in step if the map is restyled. Replacing a scale emits a ggplot2 message,
# which is expected and suppressed here.
map_rich_sankey <- suppressMessages(
  map_rich + scale_fill_manual(values = map_colors_sankey, name = "Link category")
)

map_rich_sankey

# ggsave("results/figures/map_richest_site_sankey_colors.png",
#        plot = map_rich_sankey, width = 14, height = 7, dpi = 300, bg = "white")
#
# Unlabeled counterpart, matching map_rich_blank:
# map_rich_sankey_blank <- suppressMessages(
#   map_rich_blank + scale_fill_manual(values = map_colors_sankey, name = "Link category")
# )

# -- Richest rpi (camera) site specifically, corroborated by observations --
# Same logic as above but pinned to method == "rpi" so the richest camera
# site is shown regardless of whether obs or rpi wins the overall comparison.
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
      link_category == "weakly-supported"      ~ "Weakly-supported",
      link_category == "possibly_forbidden" & add_obs  ~ "Possibly forbidden — have evidence",
      link_category == "possibly_forbidden" & !add_obs ~ "Possibly forbidden — no evidence",
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

# ---- 10. Summary ----
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
#     - phantom          — never recorded at ANY site under this method
#                          (the interaction may be genuinely impossible, or simply missed)
#     - possibly_missing — recorded at ≥1 other site under this method,
#                          but absent here (likely a detection gap)
#   TN (true negative)   — unobserved at this site and the model correctly predicted absent.
#     - possibly_forbidden — never recorded at ANY site under this method
#                          (the model concurs: the interaction is probably impossible)
#     - locally_absent   — recorded at ≥1 other site under this method,
#                          but model correctly predicts absence here (locally unsuitable)
#
# CROSS-METHOD EVIDENCE ("have evidence"): the interaction_id appears as
#   ground_truth == 1 in the orthogonal sampling method (cameras for obs, or
#   vice versa).  "% corroborated" = have evidence / category total.

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
