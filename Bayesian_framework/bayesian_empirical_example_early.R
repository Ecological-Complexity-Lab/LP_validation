# Bayesian categorisation of predicted links: Cabrera worked example
#
# Turns the deterministic eight-category assignment into a posterior over the
# eight categories, for every candidate link at every site.
#
# Inputs : serra_marin_loo_prediction_results.csv
# Outputs: posterior table (one row per link, eight posterior columns)
#          confidence maps (richest site, all sites, single-category gradient)
#          phase space over the two free rates
#          table of camera corroboration by category

library(tidyverse)

# ---- 1. Setup and parameters ----

## Figures are written to Bayesian_framework/bayesian_figures. Resolved here so
## the script works whether it is run from the project root or from inside
## Bayesian_framework/, and the folder is created if it is missing.
fig_dir <- if (dir.exists("Bayesian_framework"))
             file.path("Bayesian_framework", "bayesian_figures") else "bayesian_figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

DATA <- "results/predictions/serra_marin_loo_prediction_results.csv"

# --- Free parameters. Change these to explore. ---
RHO   <- 0.40   # realisation rate: how often a realisable link is actually
                # realised at a site. MADE UP for now. This is the single
                # least-constrained number in the analysis.
F_POS <- 0.05   # false-detection rate. Fixed a priori because nothing in the
                # data identifies it.

# --- Prior knobs. The prior for each bit runs linearly from MIN (two
#     specialists) to MAX (two generalists), driven by regional degree. ---
PI_Y <- c(min = 0.10, max = 0.70)   # would a truth-fitted model predict it?
PI_L <- c(min = 0.05, max = 0.50)   # is it realised here?
PI_R <- c(min = 0.10, max = 0.80)   # is it realisable in the replicates?

# Set to TRUE for a uniform prior (1/8 each) as the reference run.
UNIFORM_PRIOR <- FALSE

# The eight categories and their error-free evidence signatures (Table S2).
# zY = model prediction, zl = realised locally, zr = realisable in replicates.
cats <- tibble(
  category = c("recurrent", "locally unique", "possibly missing", "phantom",
               "model-elusive", "weakly supported", "locally absent",
               "possibly forbidden"),
  zY = c(1, 1, 1, 1, 0, 0, 0, 0),
  zl = c(1, 1, 0, 0, 1, 1, 0, 0),
  zr = c(1, 0, 1, 0, 1, 0, 1, 0)
)

# ---- 2. Build the evidence ----
# Each link needs four things: Y (predicted?), O_l (seen here?), n (how many
# replicate sites recorded it), R (how many replicate sites could have).
#
# Note on R. We treat every other site as a replicate, and let the realisation
# rate carry co-occurrence: rho_ij(k) = c_ij(k) * r_ij, where c_ij(k) is 1 only
# where both partners are present. At a site where they do not co-occur, both
# p1 and p0 are zero, so that site contributes an identical factor of 1 to
# every category and cancels. The arithmetic is therefore the same as counting
# only co-occurring sites, which is what R below holds. This matters: 582 of
# 1624 links co-occur nowhere else, so their R is 0 and the replicate evidence
# is silent rather than negative.

raw <- read_csv(DATA, show_col_types = FALSE)

obs <- raw %>%
  filter(method == "obs") %>%
  mutate(pair = paste(higher_level, lower_level, sep = "||"))

evidence <- obs %>%
  group_by(pair) %>%
  mutate(
    n_total  = sum(ground_truth),          # sites recording it, incl. focal
    R_total  = n_distinct(focal_site)      # sites where the pair co-occurs
  ) %>%
  ungroup() %>%
  transmute(
    link_ID, focal_site, pair,
    pollinator = higher_level,
    plant      = lower_level,
    q          = probability,              # softImpute score, kept for later
    Y          = prediction,               # model prediction
    O_l        = ground_truth,             # local observation
    n          = n_total - ground_truth,   # detections in the replicates
    R          = R_total - 1,              # replicates that could detect it
    O_r        = as.integer(n > 0)
  ) %>%
  # deterministic label, for comparison with the posterior
  left_join(cats, by = character()) %>%
  filter(Y == zY, O_l == zl, O_r == zr) %>%
  select(-zY, -zl, -zr) %>%
  rename(det_category = category)

# ---- 3. Error rates ----

# --- eps_Y: model error, proxied by 1 - F1. -------------------------------
# We use F1 rather than 1 - accuracy because the dataset has no verified true
# negatives: an unobserved pair is unknown, not absent. Accuracy and
# specificity both credit the model for 1077 unverified zeros. F1 ignores
# them, which is the right call under positive-unlabelled data.
# Caveat to state in the text: F1 is not a probability, so 1 - F1 is a
# convention, not an estimate. It is symmetric (one rate, both directions).

tp <- sum(obs$ground_truth == 1 & obs$prediction == 1)
fp <- sum(obs$ground_truth == 0 & obs$prediction == 1)
fn <- sum(obs$ground_truth == 1 & obs$prediction == 0)

precision <- tp / (tp + fp)
recall    <- tp / (tp + fn)
F1        <- 2 * precision * recall / (precision + recall)
EPS_Y     <- 1 - F1

# --- eps_l: local miss rate, from the camera contrast. ---------------------
# On site-pairs surveyed by both methods, take the camera's records as the
# reference and ask how often direct observation missed them. This is a lower
# bound, because the cameras miss links too.

cam <- raw %>%
  filter(method == "rpi") %>%
  transmute(higher_level, lower_level, focal_site, cam = ground_truth)

overlap <- obs %>%
  inner_join(cam, by = c("higher_level", "lower_level", "focal_site"))

EPS_L <- 1 - with(overlap, sum(cam == 1 & ground_truth == 1) / sum(cam == 1))

# --- derived per-replicate detection probabilities (Eq. S6) ----------------
p1 <- RHO * (1 - EPS_L)   # realisable link: realised AND detected
p0 <- F_POS               # unrealisable link: false detection only

cat(sprintf("eps_Y = %.3f  (F1 = %.3f)\n", EPS_Y, F1))
cat(sprintf("eps_l = %.3f  (from %d shared site-pairs)\n", EPS_L, nrow(overlap)))
cat(sprintf("f     = %.3f  (fixed)\n", F_POS))
cat(sprintf("rho   = %.3f  (made up)   ->  p1 = %.3f, p0 = %.3f\n",
            RHO, p1, p0))

# ---- 4. The prior, from regional degree ----
# Degree is read off the network, so it is already part of the evidence. To
# keep prior and evidence disjoint we compute each species' degree from the
# five NON-focal sites only. Same blocking logic as the leave-one-out
# prediction, so it costs nothing.
#
# For each link we form a generality score g in [0,1] from the two regional
# degrees, then read each of the three prior bits off a line from MIN to MAX.
# The three bits are assumed independent, so the prior over the eight
# categories is the product of the three.

deg_by_site <- obs %>%
  filter(ground_truth == 1) %>%
  distinct(focal_site, higher_level, lower_level)

# regional degree of each species EXCLUDING the focal site
regional_degree <- evidence %>%
  distinct(focal_site, pollinator, plant) %>%
  rowwise() %>%
  mutate(
    d_pol = deg_by_site %>%
      filter(focal_site != .env$focal_site, higher_level == .env$pollinator) %>%
      pull(lower_level) %>% n_distinct(),
    d_pla = deg_by_site %>%
      filter(focal_site != .env$focal_site, lower_level == .env$plant) %>%
      pull(higher_level) %>% n_distinct()
  ) %>%
  ungroup()

prior_tbl <- regional_degree %>%
  mutate(
    # scale each degree to [0,1], then multiply: g is high only when BOTH
    # partners are generalists
    g    = (d_pol / max(d_pol)) * (d_pla / max(d_pla)),
    pi_Y = PI_Y["min"] + (PI_Y["max"] - PI_Y["min"]) * g,
    pi_l = PI_L["min"] + (PI_L["max"] - PI_L["min"]) * g,
    pi_r = PI_R["min"] + (PI_R["max"] - PI_R["min"]) * g
  )

evidence <- evidence %>%
  left_join(prior_tbl, by = c("focal_site", "pollinator", "plant"))

# ---- 5. The posterior ----
# For each category c we ask: if c were true, how likely is the evidence we
# saw? Each source contributes (1 - eps) when it agrees with the category's
# signature and eps when it disagrees. The three contributions multiply.
# Multiply by the prior, then normalise across the eight so they sum to 1.

posterior_for <- function(zY, zl, zr, ev) {

  # model: agrees with the signature or not
  f_Y <- ifelse(ev$Y == zY, 1 - EPS_Y, EPS_Y)

  # local observation: two directions of error, a miss (eps_l) and a false
  # detection (f). Which one applies depends on the category, not the evidence.
  f_l <- case_when(
    zl == 1 & ev$O_l == 1 ~ 1 - EPS_L,   # realised and seen
    zl == 1 & ev$O_l == 0 ~ EPS_L,       # realised but missed
    zl == 0 & ev$O_l == 1 ~ F_POS,       # not realised but recorded
    zl == 0 & ev$O_l == 0 ~ 1 - F_POS    # not realised and not recorded
  )

  # replicates: n successes out of R, with a per-replicate probability that
  # depends only on whether the category says the link is realisable there.
  # The binomial coefficient is identical for all eight and cancels, so n is a
  # sufficient summary. When R = 0 this whole term is 1 and drops out.
  p   <- ifelse(zr == 1, p1, p0)
  f_r <- p^ev$n * (1 - p)^(ev$R - ev$n)

  # prior for this category, from the three independent bits
  pri <- if (UNIFORM_PRIOR) 1/8 else {
    (ev$pi_Y^zY * (1 - ev$pi_Y)^(1 - zY)) *
    (ev$pi_l^zl * (1 - ev$pi_l)^(1 - zl)) *
    (ev$pi_r^zr * (1 - ev$pi_r)^(1 - zr))
  }

  f_Y * f_l * f_r * pri
}

# unnormalised posterior, one column per category
unnorm <- map(seq_len(nrow(cats)), function(i) {
  posterior_for(cats$zY[i], cats$zl[i], cats$zr[i], evidence)
}) %>%
  set_names(cats$category) %>%
  as_tibble()

posterior <- unnorm / rowSums(unnorm)

result <- bind_cols(evidence, posterior) %>%
  mutate(
    # aggregate quantities from the main text
    phi = 1 - (phantom + `possibly forbidden`),      # feasibility confidence
    top_category = cats$category[max.col(as.matrix(posterior))],
    top_conf     = apply(as.matrix(posterior), 1, max),
    # how ambiguous is this link? 0 = certain, 1 = completely undecided
    entropy = apply(as.matrix(posterior), 1,
                    function(p) -sum(p[p > 0] * log(p[p > 0])) / log(8))
  )

write_csv(result, "bayesian_categorisation_results.csv")

# ---- 6. Summaries and figures ----

# -- 6a. Deterministic counts vs expected counts --
# The deterministic version assigns each link to exactly one box. The Bayesian
# version spreads it. Summing the posterior column gives the EXPECTED number
# of links in each category, which is the honest headline number.

comparison <- result %>%
  count(det_category, name = "deterministic") %>%
  left_join(
    posterior %>% summarise(across(everything(), sum)) %>%
      pivot_longer(everything(), names_to = "det_category",
                   values_to = "expected"),
    by = "det_category"
  ) %>%
  mutate(expected = round(expected, 1)) %>%
  arrange(desc(deterministic))

print(comparison)

# -- 6b. Confidence maps --
#
# Three plots:
#   map_rich      the richest site, hue = category, ALPHA = posterior
#   map_rich_grad the richest site, single category, white -> colour
#   map_all       all six sites, each showing only its own species
#
# map_rich is the one worth building the figure around. It keeps the exact
# palette, species ordering and layout of the published categorical map, so a
# reader who knows that figure reads this one instantly. The only change is
# that each cell now fades in proportion to how confident we are in the label
# it carries. Solid cells are well-evidenced. Ghostly cells are guesses wearing
# a category name.

# -- Palette, from serra_marin_link_classification_clean.R (Section 2). The
#    published map splits phantom and possibly forbidden by whether the other
#    sampling method corroborated them; here the two pairs are unified back
#    into single categories, because the posterior already carries that
#    evidence and a second visual channel for it would double-count. The
#    unified categories keep their no-evidence colours; the two aquamarine
#    have-evidence tints therefore drop out, leaving the eight categories of
#    the taxonomy. --
map_display_levels <- c(
  "Locally unique", "Recurrent",
  "Phantom",
  "Possibly missing", "Weakly-supported", "Model elusive",
  "Possibly forbidden",
  "Locally absent"
)

map_colors <- c(
  # -- Observed in >=1 replicate — warm --
  "Recurrent"          = "#C05030",
  "Model elusive"      = "#8A3050",
  "Possibly missing"   = "#DBA040",
  "Locally absent"     = "#F9C8B0",
  # -- Not observed in any replicate — cool --
  "Locally unique"     = "#B0C4DE",
  "Weakly-supported"   = "#E8A8B8",
  "Phantom"            = "#D5C8F0",
  "Possibly forbidden" = "#5E3DA0"
)

## --- 1. attach the confidence in each link's OWN category ------------------
# For every link, pull the posterior column that matches its deterministic
# label. This single number is what drives the transparency, for every one of
# the eight categories.

post_mat <- as.matrix(result[, cats$category])

map_df <- result %>%
  mutate(
    conf_det = post_mat[cbind(seq_len(nrow(post_mat)),
                              match(det_category, cats$category))],
    link_display = factor(case_when(
      det_category == "recurrent"          ~ "Recurrent",
      det_category == "locally unique"     ~ "Locally unique",
      det_category == "possibly missing"   ~ "Possibly missing",
      det_category == "phantom"            ~ "Phantom",
      det_category == "model-elusive"      ~ "Model elusive",
      det_category == "weakly supported"   ~ "Weakly-supported",
      det_category == "locally absent"     ~ "Locally absent",
      det_category == "possibly forbidden" ~ "Possibly forbidden"
    ), levels = map_display_levels)
  )

# Posterior confidence differs enormously in LEVEL between categories: possibly
# forbidden sits around 0.8 while locally unique sits around 0.1, so a single
# alpha scale renders one category solid and the rest almost invisible, and the
# hues stop being readable. Rescaling within each category puts every category
# on its own 0-1 span, so alpha reads as "confident FOR THIS CATEGORY" and the
# colours stay legible. The cost is that alpha is no longer comparable across
# categories, which the legend has to say. Categories whose cells all share one
# value collapse to the top of the range rather than dividing by zero.
rescale_within <- function(x) {
  r <- range(x)
  if (diff(r) < 1e-12) rep(1, length(x)) else (x - r[1]) / diff(r)
}

## --- 2. species ordering ---------------------------------------------------
# Matched to map_rich in serra_marin_link_classification_clean.R (Section 9).
# Degree is the number of DISTINCT partners recorded under this method, pooled
# across all sites, so a generalist sits in the same place in every panel. Note
# this is not the same statistic as counting (site, partner) rows, which sums
# the per-site degrees and gives a different order.
overall_poll_degree <- obs %>%
  filter(ground_truth == 1) %>%
  group_by(higher_level) %>%
  summarise(overall_poll_degree = n_distinct(lower_level), .groups = "drop")

overall_plant_degree <- obs %>%
  filter(ground_truth == 1) %>%
  group_by(lower_level) %>%
  summarise(overall_plant_degree = n_distinct(higher_level), .groups = "drop")

# Levels descending in degree, and NOT reversed on the y axis: ggplot places
# the first level at the bottom, which is what the published map does.
pol_order <- overall_poll_degree  %>% arrange(desc(overall_poll_degree))  %>% pull(higher_level)
pla_order <- overall_plant_degree %>% arrange(desc(overall_plant_degree)) %>% pull(lower_level)

map_df <- map_df %>%
  mutate(
    pollinator = factor(pollinator, levels = pol_order),
    plant      = factor(plant,      levels = pla_order)
  )

# italic species labels, as in the published map
ital <- function(x) lapply(strsplit(x, " "),
                           function(y) bquote(italic(.(paste(y, collapse = " ")))))

## --- 3. the richest site ---------------------------------------------------

richness <- map_df %>%
  filter(O_l == 1) %>%
  group_by(focal_site) %>%
  summarise(n_plants      = n_distinct(plant),
            n_pollinators = n_distinct(pollinator),
            n_species     = n_plants + n_pollinators,
            .groups = "drop") %>%
  arrange(desc(n_species))

print(richness)

RICH_SITE <- richness$focal_site[1]     # "Rocky coastal 1": 11 plants, 46 pollinators
cat(sprintf("\nRichest site: %s\n", RICH_SITE))

# droplevels() is what restricts the axes to species present at this site.
# The candidate set is the full plant x pollinator grid per site, so the panel
# comes out as a complete rectangle with no empty rows or columns.
rich_df <- map_df %>% filter(focal_site == RICH_SITE) %>% droplevels()

# Re-level from the site's own rows, which is what the reference map_rich does.
# Several species tie on degree, and a tie is broken by the order the rows come
# in: taking distinct() over the whole dataset and then subsetting gives a
# different winner than taking it over this site alone. Ordering here reproduces
# the published axis exactly.
pol_order_rich <- rich_df %>%
  distinct(pollinator) %>%
  left_join(overall_poll_degree, by = c("pollinator" = "higher_level")) %>%
  arrange(desc(overall_poll_degree)) %>% pull(pollinator) %>% as.character()

pla_order_rich <- rich_df %>%
  distinct(plant) %>%
  left_join(overall_plant_degree, by = c("plant" = "lower_level")) %>%
  arrange(desc(overall_plant_degree)) %>% pull(plant) %>% as.character()

rich_df <- rich_df %>%
  mutate(pollinator = factor(as.character(pollinator), levels = pol_order_rich),
         plant      = factor(as.character(plant),      levels = pla_order_rich))

## --- 4. MAP_RICH: category hue, posterior transparency ---------------------
# Alpha is the posterior in the cell's own category, rescaled within that
# category (see rescale_within above). The floor is well above zero so that the
# least confident cell still shows its hue: the point of the figure is to read
# category and confidence together, which fails if half the map is white.

ALPHA_MAX <- max(rich_df$conf_det)
cat(sprintf("raw posterior at this site spans %.2f to %.2f; alpha is rescaled within category\n",
            min(rich_df$conf_det), ALPHA_MAX))

rich_df <- rich_df %>%
  group_by(link_display) %>%
  mutate(conf_rel = rescale_within(conf_det)) %>%
  ungroup()

map_rich <- ggplot(rich_df, aes(x = pollinator, y = plant,
                                fill = link_display, alpha = conf_rel)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = map_colors, name = "Link category", drop = FALSE) +
  scale_alpha_continuous(range = c(0.35, 1), limits = c(0, 1),
                         breaks = c(0, 0.5, 1),
                         labels = c("lowest", "mid", "highest"),
                         name = "Posterior\nconfidence\n(within category)") +
  # keep the category legend on top, as in the published map
  guides(fill  = guide_legend(override.aes = list(alpha = 1), order = 1),
         alpha = guide_legend(order = 2)) +
  labs(x = "Pollinator", y = "Plant") +
  scale_x_discrete(labels = ital) +
  scale_y_discrete(labels = ital) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 13, angle = 90, vjust = 0.5),
    axis.text.y  = element_text(size = 13),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    panel.grid   = element_blank()
  )

map_rich
ggsave(file.path(fig_dir, "fig_map_rich_confidence.pdf"), map_rich,
       width = 14, height = 7)

## --- 5. MAP_RICH_GRAD: one category, white to colour ----------------------
# Shai's original request. Cleaner for a single message, and the only version
# that shows the posterior for a category a link was NOT assigned to.

FOCAL_CATEGORY <- "possibly missing"

map_rich_grad <- rich_df %>%
  mutate(conf = .data[[FOCAL_CATEGORY]]) %>%
  ggplot(aes(x = pollinator, y = plant, fill = conf)) +
  geom_tile(colour = "grey90", linewidth = 0.4) +
  scale_fill_gradient(low = "white", high = unname(map_colors["Possibly missing"]),
                      limits = c(0, ALPHA_MAX), name = "Posterior") +
  labs(x = "Pollinator", y = "Plant",
       title = paste0(RICH_SITE, ": confidence that each link is ",
                      FOCAL_CATEGORY)) +
  scale_x_discrete(labels = ital) +
  scale_y_discrete(labels = ital) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 11),
    panel.grid  = element_blank()
  )

map_rich_grad
ggsave(file.path(fig_dir, "fig_map_rich_gradient.pdf"), map_rich_grad, width = 14, height = 7)

## --- 6. MAP_ALL: all six sites, each with its own species -----------------
# scales = "free" is what lets each panel drop the species it does not have.
# Panels come out at different sizes; space = "free" in facet_grid would fix
# that but only for a single row, so facet_wrap with free scales is the
# practical choice for six panels.

# Rescaled within category across all six sites together, so a category means
# the same shade in every panel and the panels stay comparable.
map_df <- map_df %>%
  group_by(link_display) %>%
  mutate(conf_rel = rescale_within(conf_det)) %>%
  ungroup()

map_all <- ggplot(map_df, aes(x = pollinator, y = plant,
                              fill = link_display, alpha = conf_rel)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  facet_wrap(~ focal_site, ncol = 2, scales = "free") +
  scale_fill_manual(values = map_colors, name = "Link category", drop = FALSE) +
  scale_alpha_continuous(range = c(0.35, 1), limits = c(0, 1),
                         breaks = c(0, 0.5, 1),
                         labels = c("lowest", "mid", "highest"),
                         name = "Posterior\nconfidence\n(within category)") +
  guides(fill  = guide_legend(override.aes = list(alpha = 1), order = 1),
         alpha = guide_legend(order = 2)) +
  labs(x = "Pollinator", y = "Plant") +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 4, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 5),
    panel.grid  = element_blank()
  )

ggsave(file.path(fig_dir, "fig_map_all_sites.pdf"), map_all, width = 16, height = 20, limitsize = FALSE)

# --- 6c. FIGURE 2: phase space over the two free rates --------------------
# Shai's second suggestion. eps_l and f are the two rates we are least sure
# of, so sweep them and watch a single summary move. Here: how many links
# reach TAU confidence of being possibly missing.
#
# IMPORTANT: TAU must sit below kappa, the maximum contextual confidence. With
# eps_Y = 1 - F1 = 0.417 and eps_l = 0.34, kappa is only about 0.43, so a
# threshold of 0.5 would return zero links everywhere and the panel would be
# blank. Check kappa before choosing TAU.

TAU <- 0.35
kappa <- (1 - EPS_Y) * (1 - F_POS) / ((1 - F_POS) + EPS_L)
cat(sprintf("kappa (ceiling when O_l = 0) = %.3f ; TAU = %.2f\n", kappa, TAU))

sweep <- expand_grid(el = seq(0.05, 0.60, by = 0.025),
                     ff = seq(0.01, 0.20, by = 0.01)) %>%
  rowwise() %>%
  mutate(n_links = {
    EPS_L <<- el; F_POS <<- ff
    p1 <<- RHO * (1 - el); p0 <<- ff
    u <- map(seq_len(nrow(cats)), function(i)
      posterior_for(cats$zY[i], cats$zl[i], cats$zr[i], evidence)) %>%
      set_names(cats$category) %>% as_tibble()
    sum((u / rowSums(u))[["possibly missing"]] > TAU)
  }) %>%
  ungroup()

fig2 <- ggplot(sweep, aes(el, ff, fill = n_links)) +
  geom_raster() +
  scale_fill_viridis_c(name = paste0("Links > ", TAU)) +
  labs(x = expression(epsilon[l]~"(local miss rate)"),
       y = "f (false detection rate)",
       title = paste0("Links reaching ", TAU, " confidence of being possibly missing")) +
  theme_minimal()

ggsave(file.path(fig_dir, "fig_phase_space.pdf"), fig2, width = 7, height = 5)

# restore the chosen rates after the sweep
EPS_L <- 1 - with(overlap, sum(cam == 1 & ground_truth == 1) / sum(cam == 1))
F_POS <- 0.05
p1 <- RHO * (1 - EPS_L); p0 <- F_POS

# --- 6d. Camera corroboration by category ---------------------------------
# The independent check on the categorisation. For links NOT seen locally by
# direct observation, ask how often the cameras found them, split by category.
#
# This is deliberately deterministic. It uses no error rates and no posterior,
# so it tests the categorisation itself rather than the assumed inputs. We do
# not calibrate the posterior against the cameras: the cameras already set
# eps_l, so grading the posterior with them would be circular, and only 34 of
# the 51 possibly missing links have camera coverage anyway.

corroboration <- result %>%
  filter(O_l == 0) %>%                       # links direct observation missed
  inner_join(cam, by = c("pollinator" = "higher_level",
                         "plant"      = "lower_level",
                         "focal_site" = "focal_site")) %>%
  group_by(det_category) %>%
  summarise(
    n_links   = n(),
    confirmed = sum(cam == 1),
    rate      = mean(cam == 1),
    .groups   = "drop"
  ) %>%
  arrange(desc(rate))

print(corroboration)
write_csv(corroboration, "camera_corroboration_by_category.csv")

cat("\nDone. Wrote two CSVs and two figures.\n")
