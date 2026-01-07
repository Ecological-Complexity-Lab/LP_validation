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

## ---- no link withholding ----
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
    is_shared     = n_obs_total >= 2,  # with 2 methods, this means == 2
    obs_elsewhere = (n_obs_total - (original_binary == 1)) >= 1,
    
    validation = ifelse(obs_elsewhere, "Corroborated", "Not corroborated")
  )

# classify links
df_categorized <- df_flagged %>%
  mutate(
    link_category = case_when(
      # observed in only this method
      is_unique   & original_binary == 1 & predicted_bin == 1 ~ "locally_unique_links", # method specific links
      is_unique   & original_binary == 1 & predicted_bin == 0 ~ "unsupported_links",
      
      # observed in both methods
      is_shared   & original_binary == 1 & predicted_bin == 1 ~ "confirmed_links", # recurrent
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

# add evaluation classes
df_categorized <- df_categorized %>%
  mutate(
    confusion = case_when(
      original_binary == 1 & predicted_bin == 1 ~ "TP",
      original_binary == 1 & predicted_bin == 0 ~ "FN",
      original_binary == 0 & predicted_bin == 1 ~ "FP",
      original_binary == 0 & predicted_bin == 0 ~ "TN",
      TRUE ~ "UNK"
    )
  )

# summary table for plot
summary_tbl <- df_categorized %>%
  filter(link_category != "unclassified", confusion != "UNK") %>%
  count(method, confusion, validation, link_category, name = "value")

## ---- alluvial plot ----
# make an alluvial plot
method_to_plot <- "method2"  # we can change to "method2"

summary_tbl1 <- summary_tbl %>%
  filter(method == method_to_plot)

cat_means <- summary_tbl %>%
  filter(method == method_to_plot) %>%
  group_by(link_category) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  tibble::deframe()   # named vector: names are link_category

# a helper for the plot
get_mean <- function(nm) if (nm %in% names(cat_means)) unname(cat_means[[nm]]) else 0

# set flows
flows_tbl <- summary_tbl1 %>%
  mutate(
    L1 = confusion,
    L2 = validation,
    L3 = link_category
  ) %>%
  select(method, L1, L2, L3, value) %>%
  group_by(method, L1, L2, L3) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  group_by(method) %>%
  mutate(
    total_confusion = sum(value),
    prop = ifelse(total_confusion > 0, value / total_confusion, 0),
    alluvium_id = paste(L1, L3, sep = "⟂")
  ) %>%
  ungroup()

# Aesthetics — tweak these as you like
col_confusion <- c(
  "TP"="lightsteelblue","FP"="lightsteelblue2","TN"="rosybrown","FN"="rosybrown2"
)
col_validation <- c(
  "Corroborated"="sandybrown","Not corroborated"="thistle3"
)
col_subcats <- c(
  "confirmed_links"="coral3","locally_unique_links"="thistle3",
  "cryptic_links"="coral","unsupported_links"="thistle1",
  "possibly_missing_links"="coral2","spurious_links"="thistle",
  "feasible_links"="coral1","likely_forbidden"="thistle2"
)
stratum_fill <- c(col_confusion, col_validation, col_subcats)

flow_alpha <- 0.7     # flow transparency
flow_colour <- NA      # outline color for flows; e.g., "grey30" to draw borders
stratum_label_size <- 4

# Choose whether to plot mean counts or proportions:
metric <- "prop"  # set to "value" for raw mean counts

# Prepare "lodes" format for 3 axes and KEEP alluvium_id
flows_long <- flows_tbl %>%
  select(L1, L2, L3, value, prop, alluvium_id) %>%
  pivot_longer(cols = c(L1, L2, L3),
               names_to = "axis", values_to = "stratum") %>%
  mutate(
    axis = recode(axis, L1="Confusion", L2="Validation", L3="Subtype"),
    axis = factor(axis, levels = c("Confusion","Validation","Subtype"))
  )

# Pretty labels for strata
stratum_labeller <- function(x) {
  x %>% str_replace_all("_", " ") %>% str_to_sentence()
}

## -------- choose your orders here --------
order_confusion  <- c("TP", "FP", "TN", "FN")                 # LEFT column order
order_validation <- c("Not corroborated", "Corroborated") # MIDDLE column order
order_subtypes   <- c(                                       # RIGHT column order
  "locally_unique_links", "spurious_links",
  "likely_forbidden", "unsupported_links",
  "confirmed_links", "possibly_missing_links",
  "feasible_links", "cryptic_links"
)

# Rebuild flows_long with your custom orders applied
flows_long <- flows_tbl %>%
  select(L1, L2, L3, value, prop, alluvium_id) %>%
  tidyr::pivot_longer(c(L1, L2, L3), names_to = "axis", values_to = "stratum") %>%
  dplyr::mutate(
    axis = dplyr::recode(axis, L1 = "Confusion", L2 = "Validation", L3 = "Subtype"),
    # make stratum a factor with the order you chose, depending on axis
    stratum = dplyr::case_when(
      axis == "Confusion"  ~ factor(stratum, levels = order_confusion),
      axis == "Validation" ~ factor(stratum, levels = order_validation),
      axis == "Subtype"    ~ factor(stratum, levels = order_subtypes),
      TRUE ~ factor(stratum)
    ),
    # also lock the axis order (left -> middle -> right)
    axis = factor(axis, levels = c("Confusion", "Validation", "Subtype"))
  )

gg <- ggplot(
  flows_long,
  aes(x = axis,
      stratum = stratum,
      alluvium = alluvium_id,
      y = .data[[metric]],
      fill = stratum)
) +
  geom_stratum(width = 0.03, color = "white") +
  geom_alluvium(color = flow_colour, alpha = flow_alpha, width = 0.25) +
  scale_fill_manual(values = stratum_fill, guide = "none") +
  scale_y_continuous(labels = if (metric=="prop") percent_format(accuracy = 1) else label_number_si()) +
  labs(
    title = if (metric=="prop") "Alluvial of mean proportions across iterations"
    else "Alluvial of mean counts across iterations",
    subtitle = "Confusion classes → Validation group → Subcategories",
    x = NULL, y = if (metric=="prop") "Proportion of interactions" else "Mean count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

gg

bg_col <- "white"   # or "#F7F4EF" or whatever your panel background is

gg <- ggplot(
  flows_long,
  aes(x = axis,
      stratum = stratum,
      alluvium = alluvium_id,
      y = .data[[metric]],
      fill = stratum)
) +
  # 2) Wide "mask" strata: same width as flows, fill = background,
  #    so they trim the flows exactly at the axis
  geom_stratum(
    width  = 0.02,
    color  = NA,
    fill   = bg_col,
    alpha  = 1
  ) +
  
  # 3) Narrow visible strata on top
  geom_stratum(
    width = 0.03,
    color = "white"
  ) +
  # 1) Flows first
  geom_alluvium(
    color = flow_colour,
    alpha = flow_alpha,
    width = 0.15,
    knot.pos = 0.2
  ) +
  
  scale_fill_manual(values = stratum_fill, guide = "none") +
  scale_y_continuous(
    labels = if (metric == "prop") percent_format(accuracy = 1)
    else label_number_si()
  ) +
  labs(
    title    = if (metric=="prop") "Alluvial of mean proportions across iterations"
    else "Alluvial of mean counts across iterations",
    subtitle = "Confusion classes → Validation group → Subcategories",
    x = NULL,
    y = if (metric=="prop") "Proportion of interactions" else "Mean count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 12, face = "bold"),
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(color = "grey30")
  )


# totals per stratum (for % text)
stratum_totals <- flows_long %>%
  dplyr::group_by(axis, stratum) %>%
  dplyr::summarise(val = sum(.data[[metric]]), .groups = "drop") %>%
  dplyr::mutate(stratum_chr = as.character(stratum))

# build one hidden stratum layer to grab box geometry (same width as above: 0.25)
tmp_build <- ggplot_build(
  ggplot(flows_long,
         aes(x = axis, stratum = stratum, alluvium = alluvium_id, y = .data[[metric]])) +
    geom_stratum(width = 0.03, color = NA)
)

geo <- as.data.frame(tmp_build$data[[1]])
ax_levels <- levels(flows_long$axis)

# midpoints and axis name for each stratum box
label_geom <- geo %>%
  dplyr::transmute(
    x_mid       = (xmin + xmax)/2,
    y_mid       = (ymin + ymax)/2,
    stratum_chr = as.character(stratum),
    axis        = ax_levels[pmax(1, pmin(length(ax_levels), round((xmin + xmax)/2)))]
  )

# join values + geometry, build the two-line label strings
label_df <- dplyr::left_join(
  stratum_totals,
  label_geom,
  by = c("axis","stratum_chr")
) %>%
  dplyr::mutate(
    name_txt = stratum_chr %>% stringr::str_replace_all("_"," ") %>% stringr::str_to_sentence(),
    pct_txt  = if (metric == "prop") scales::percent(val, accuracy = 1)
    else scales::label_number_si()(val)
  )

# tweakable label settings
label_nudge_x   <- 0.03                                 # push labels to the right of each box
label_nudge_y   <- 0.05  
y_off           <- 0.5 * diff(range(flows_long[[metric]], na.rm = TRUE))  # vertical gap for % line
name_size       <- 4.2
pct_size        <- 3.6
font_family     <- ""                                     # "" = default device font
label_color_map <- c("Confusion"="mistyrose4","Validation"="mistyrose4","Subtype"="mistyrose4")

label_df <- label_df %>%
  mutate(
    label_final = paste0(name_txt, " (", pct_txt, ")")
  )

gg <- gg +
  geom_text(
    data = label_df,
    inherit.aes = FALSE,
    aes(x = x_mid + label_nudge_x, y = y_mid, label = label_final, color = axis),
    fontface = "bold",
    size = name_size,
    family = font_family,
    hjust = 0
  ) +
  scale_color_manual(values = label_color_map, guide = "none") +
  coord_cartesian(clip = "off")

gg <- ggplot(
  flows_long,
  aes(x = axis,
      stratum = stratum,
      alluvium = alluvium_id,
      y = .data[[metric]],
      fill = stratum)
) +
  # 2) Wide "mask" strata: same width as flows, fill = background,
  #    so they trim the flows exactly at the axis
  geom_stratum(
    width  = 0.02,
    color  = NA,
    fill   = bg_col,
    alpha  = 1
  ) +
  
  # 3) Narrow visible strata on top
  geom_stratum(
    width = 0.03,
    color = "white"
  ) +
  
  # 1) Flows
  geom_alluvium(
    color    = flow_colour,
    alpha    = flow_alpha,
    width    = 0.15,
    knot.pos = 0.2
  ) +
  
  scale_fill_manual(values = stratum_fill, guide = "none") +
  scale_y_continuous(
    labels = if (metric == "prop") percent_format(accuracy = 1)
    else label_number_si()
  ) +
  labs(
    title    = NULL,   # remove title
    subtitle = NULL,   # remove subtitle
    x = NULL,
    y = if (metric == "prop") "Proportion of interactions" else "Mean count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 12, face = "bold"),
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(color = "grey30")
  )

# range for vertical offsets
yr <- diff(range(flows_long[[metric]], na.rm = TRUE))

# total value per stratum per axis
stratum_totals <- flows_long %>%
  dplyr::group_by(axis, stratum) %>%
  dplyr::summarise(val = sum(.data[[metric]]), .groups = "drop") %>%
  dplyr::mutate(stratum_chr = as.character(stratum))

# build a hidden stratum layer to get box geometry
tmp_build <- ggplot_build(
  ggplot(flows_long,
         aes(x = axis,
             stratum = stratum,
             alluvium = alluvium_id,
             y = .data[[metric]])) +
    geom_stratum(width = 0.03, color = NA)
)

geo <- as.data.frame(tmp_build$data[[1]])
ax_levels <- levels(flows_long$axis)

# midpoints and axis name for each stratum box
label_geom <- geo %>%
  dplyr::transmute(
    x_mid       = (xmin + xmax)/2,
    y_mid       = (ymin + ymax)/2,
    stratum_chr = as.character(stratum),
    axis        = ax_levels[round(x)]  # map numeric x back to axis factor
  )

# join values + geometry, build label text and vertically stagger labels
label_df <- dplyr::left_join(
  stratum_totals,
  label_geom,
  by = c("axis", "stratum_chr")
) %>%
  dplyr::mutate(
    name_txt = stratum_chr %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_to_sentence(),
    pct_txt  = if (metric == "prop") scales::percent(val, accuracy = 1)
    else scales::label_number_si()(val),
    label_final = paste0(name_txt, " (", pct_txt, ")")
  ) %>%
  dplyr::arrange(axis, y_mid) %>%
  dplyr::group_by(axis) %>%
  dplyr::mutate(
    # stagger labels within each axis to reduce overlap
    label_y = y_mid + (row_number() - mean(row_number())) * (0.033 * yr)
  ) %>%
  dplyr::ungroup()

label_nudge_x <- 0.03   # adjust if you want labels further right
label_nudge_y <- 0.02

gg <- gg +
  geom_text(
    data = label_df,
    inherit.aes = FALSE,
    aes(x = x_mid + label_nudge_x,
        y = label_y + label_nudge_y,
        label = label_final,
        color = axis),
    fontface = "bold",
    size     = name_size,
    family   = font_family,
    hjust    = 0
  ) +
  scale_color_manual(values = label_color_map, guide = "none") +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(20, 20, 20, 20)   # room on the right for labels
  )

gg

gg <- gg +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),     # remove tick labels
    axis.ticks = element_blank(),    # remove tick marks
    axis.title = element_blank(),    # remove axis titles
    axis.line = element_blank(),     # remove axis lines
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    plot.margin = margin(20, 20, 20, 20)  # keep right margin for labels
  )

gg

# no label version:

gg <- ggplot(
  flows_long,
  aes(x = axis,
      stratum = stratum,
      alluvium = alluvium_id,
      y = .data[[metric]],
      fill = stratum)
) +
  geom_stratum(width = 0.03, color = "white") +
  geom_alluvium(color = flow_colour, alpha = flow_alpha, width = 0.25) +
  scale_fill_manual(values = stratum_fill, guide = "none") +
  scale_y_continuous(labels = if (metric=="prop") percent_format(accuracy = 1) else label_number_si()) +
  labs(
    title = if (metric=="prop") "Alluvial of mean proportions across iterations"
    else "Alluvial of mean counts across iterations",
    subtitle = "Confusion classes → Validation group → Subcategories",
    x = NULL, y = if (metric=="prop") "Proportion of interactions" else "Mean count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

bg_col <- "white"   # or "#F7F4EF" or whatever your panel background is

gg <- ggplot(
  flows_long,
  aes(x = axis,
      stratum = stratum,
      alluvium = alluvium_id,
      y = .data[[metric]],
      fill = stratum)
) +
  # 2) Wide "mask" strata: same width as flows, fill = background,
  #    so they trim the flows exactly at the axis
  geom_stratum(
    width  = 0.02,
    color  = NA,
    fill   = bg_col,
    alpha  = 1
  ) +
  
  # 3) Narrow visible strata on top
  geom_stratum(
    width = 0.03,
    color = "white"
  ) +
  # 1) Flows first
  geom_alluvium(
    color = flow_colour,
    alpha = flow_alpha,
    width = 0.15,
    knot.pos = 0.2
  ) +
  
  scale_fill_manual(values = stratum_fill, guide = "none") +
  scale_y_continuous(
    labels = if (metric == "prop") percent_format(accuracy = 1)
    else label_number_si()
  ) +
  labs(
    title    = if (metric=="prop") "Alluvial of mean proportions across iterations"
    else "Alluvial of mean counts across iterations",
    subtitle = "Confusion classes → Validation group → Subcategories",
    x = NULL,
    y = if (metric=="prop") "Proportion of interactions" else "Mean count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 12, face = "bold"),
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(color = "grey30")
  )


# totals per stratum (for % text)
stratum_totals <- flows_long %>%
  dplyr::group_by(axis, stratum) %>%
  dplyr::summarise(val = sum(.data[[metric]]), .groups = "drop") %>%
  dplyr::mutate(stratum_chr = as.character(stratum))

# build one hidden stratum layer to grab box geometry (same width as above: 0.25)
tmp_build <- ggplot_build(
  ggplot(flows_long,
         aes(x = axis, stratum = stratum, alluvium = alluvium_id, y = .data[[metric]])) +
    geom_stratum(width = 0.03, color = NA)
)

geo <- as.data.frame(tmp_build$data[[1]])
ax_levels <- levels(flows_long$axis)

# midpoints and axis name for each stratum box
label_geom <- geo %>%
  dplyr::transmute(
    x_mid       = (xmin + xmax)/2,
    y_mid       = (ymin + ymax)/2,
    stratum_chr = as.character(stratum),
    axis        = ax_levels[pmax(1, pmin(length(ax_levels), round((xmin + xmax)/2)))]
  )

# join values + geometry, build the two-line label strings
label_df <- dplyr::left_join(
  stratum_totals,
  label_geom,
  by = c("axis","stratum_chr")
) %>%
  dplyr::mutate(
    name_txt = stratum_chr %>% stringr::str_replace_all("_"," ") %>% stringr::str_to_sentence(),
    pct_txt  = if (metric == "prop") scales::percent(val, accuracy = 1)
    else scales::label_number_si()(val)
  )

# tweakable label settings
label_nudge_x   <- 0.03                                 # push labels to the right of each box
label_nudge_y   <- 0.05  
y_off           <- 0.5 * diff(range(flows_long[[metric]], na.rm = TRUE))  # vertical gap for % line
name_size       <- 4.2
pct_size        <- 3.6
font_family     <- ""                                     # "" = default device font
label_color_map <- c("Confusion"="mistyrose4","Validation"="mistyrose4","Subtype"="mistyrose4")

label_df <- label_df %>%
  mutate(
    label_final = paste0(name_txt, " (", pct_txt, ")")
  )

gg <- gg +
  # geom_text(
  #   data = label_df,
  #   inherit.aes = FALSE,
  #   aes(x = x_mid + label_nudge_x, y = y_mid, label = label_final, color = axis),
  #   fontface = "bold",
  #   size = name_size,
  #   family = font_family,
  #   hjust = 0
  # ) +
  scale_color_manual(values = label_color_map, guide = "none") +
  coord_cartesian(clip = "off")

gg <- ggplot(
  flows_long,
  aes(x = axis,
      stratum = stratum,
      alluvium = alluvium_id,
      y = .data[[metric]],
      fill = stratum)
) +
  # 2) Wide "mask" strata: same width as flows, fill = background,
  #    so they trim the flows exactly at the axis
  geom_stratum(
    width  = 0.02,
    color  = NA,
    fill   = bg_col,
    alpha  = 1
  ) +
  
  # 3) Narrow visible strata on top
  geom_stratum(
    width = 0.03,
    color = "white"
  ) +
  
  # 1) Flows
  geom_alluvium(
    color    = flow_colour,
    alpha    = flow_alpha,
    width    = 0.15,
    knot.pos = 0.2
  ) +
  
  scale_fill_manual(values = stratum_fill, guide = "none") +
  scale_y_continuous(
    labels = if (metric == "prop") percent_format(accuracy = 1)
    else label_number_si()
  ) +
  labs(
    title    = NULL,   # remove title
    subtitle = NULL,   # remove subtitle
    x = NULL,
    y = if (metric == "prop") "Proportion of interactions" else "Mean count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 12, face = "bold"),
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(color = "grey30")
  )

# range for vertical offsets
yr <- diff(range(flows_long[[metric]], na.rm = TRUE))

# total value per stratum per axis
stratum_totals <- flows_long %>%
  dplyr::group_by(axis, stratum) %>%
  dplyr::summarise(val = sum(.data[[metric]]), .groups = "drop") %>%
  dplyr::mutate(stratum_chr = as.character(stratum))

# build a hidden stratum layer to get box geometry
tmp_build <- ggplot_build(
  ggplot(flows_long,
         aes(x = axis,
             stratum = stratum,
             alluvium = alluvium_id,
             y = .data[[metric]])) +
    geom_stratum(width = 0.03, color = NA)
)

geo <- as.data.frame(tmp_build$data[[1]])
ax_levels <- levels(flows_long$axis)

# midpoints and axis name for each stratum box
label_geom <- geo %>%
  dplyr::transmute(
    x_mid       = (xmin + xmax)/2,
    y_mid       = (ymin + ymax)/2,
    stratum_chr = as.character(stratum),
    axis        = ax_levels[round(x)]  # map numeric x back to axis factor
  )

# join values + geometry, build label text and vertically stagger labels
label_df <- dplyr::left_join(
  stratum_totals,
  label_geom,
  by = c("axis", "stratum_chr")
) %>%
  dplyr::mutate(
    name_txt = stratum_chr %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_to_sentence(),
    pct_txt  = if (metric == "prop") scales::percent(val, accuracy = 1)
    else scales::label_number_si()(val),
    label_final = paste0(name_txt, " (", pct_txt, ")")
  ) %>%
  dplyr::arrange(axis, y_mid) %>%
  dplyr::group_by(axis) %>%
  dplyr::mutate(
    # stagger labels within each axis to reduce overlap
    label_y = y_mid + (row_number() - mean(row_number())) * (0.033 * yr)
  ) %>%
  dplyr::ungroup()

label_nudge_x <- 0.03   # adjust if you want labels further right
label_nudge_y <- 0.02

gg <- gg +
  # geom_text(
  #   data = label_df,
  #   inherit.aes = FALSE,
  #   aes(x = x_mid + label_nudge_x,
  #       y = label_y + label_nudge_y,
  #       label = label_final,
  #       color = axis),
  #   fontface = "bold",
  #   size     = name_size,
  #   family   = font_family,
  #   hjust    = 0
# ) +
scale_color_manual(values = label_color_map, guide = "none") +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(20, 20, 20, 20)   # room on the right for labels
  )

gg

gg <- gg +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),     # remove tick labels
    axis.ticks = element_blank(),    # remove tick marks
    axis.title = element_blank(),    # remove axis titles
    axis.line = element_blank(),     # remove axis lines
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    plot.margin = margin(20, 20, 20, 20)  # keep right margin for labels
  )

gg

# 
# # summary table
# final_table_summary <- df_categorized %>%
#   group_by(Category) %>%
#   summarise(
#     mean = mean(Count), sd = sd(Count),
#     min = min(Count), max = max(Count),
#     .groups = "drop"
#   ) %>%
#   arrange(Category)

# and a heatmap marking where are the missing links

## ---- venn ----
