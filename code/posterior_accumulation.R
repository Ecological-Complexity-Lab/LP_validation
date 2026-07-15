# ============================================================================
# Posterior certainty in the CLASS ASSIGNMENT — accumulation curves
# Case study: Serra-Marin et al. (2025), Cabrera plant-pollinator direct
#             observations, with leave-one-out link predictions.
#
# AIM (per the Box): show how confident we are that a link's assigned
# taxonomic CLASS is correct, and how that confidence accumulates as we
# fold in replicate evidence one network at a time.
#
# FORMULATION (Option A, two latent states):
#   * theta lives on just two latent states: real vs not-real.
#   * Each of the 8 classes stakes a claim on ONE of those states
#     (e.g. "possibly missing" claims the link is real; "phantom" claims
#     it is not real). See `claims_real` / `claims_notreal` below.
#   * The two-state posterior P(real | evidence) is updated sequentially.
#   * Certainty in the class assignment = the posterior mass on the state
#     the assigned class claims:
#         if class claims real     -> certainty = P(real | evidence)
#         if class claims not-real -> certainty = 1 - P(real | evidence)
#   This is a class-assignment posterior, NOT a link-existence probability.
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# ----------------------------------------------------------------------------
# PARAMETERS  (the two modelling choices; both ILLUSTRATIVE, see notes)
# ----------------------------------------------------------------------------
# THETA_REAL : P(observe the interaction in one co-occurring network | link is
#              REAL). Set well below 1 because real interactions are routinely
#              undersampled locally (detection is imperfect; Peralta et al.).
#              Not estimated here; an illustrative value. Could instead be
#              estimated from high-confidence (recurrent) links, or swept over
#              a range as a sensitivity check.
THETA_REAL <- 0.35

# EPS : P(a spurious observation in one network | link is NOT real). Small,
#       reflecting rare false detections / misidentifications. Illustrative.
EPS <- 0.03

# ----------------------------------------------------------------------------
# 1. LOAD AND RESHAPE
# ----------------------------------------------------------------------------
# CSV COLUMNS (input data, one row per candidate pair at one focal site):
#   higher_level : pollinator species
#   lower_level  : plant species
#   ground_truth : was the interaction OBSERVED at this focal site? (1/0)
#                  -> this is O_local for the focal row, and the replicate
#                     evidence when read from the OTHER co-occurring sites.
#   probability  : model's predicted probability for this link  -> Yhat (cont.)
#   prediction   : probability thresholded to 0/1               -> Yhat_bin
#   focal_site   : which of the 6 sites is the held-out focal network
#   method       : 'obs' (direct observation) or 'rpi' (camera); we use 'obs'
df <- read.csv("results/predictions/serra_martin_loo_prediction_results.csv", stringsAsFactors = FALSE)
d  <- df %>% filter(method == "obs")

# pair : unique species-pair ID, so we can collect a pair's behaviour across
#        all sites where the two species co-occur (i.e. where it is a candidate)
d <- d %>% mutate(pair = paste(higher_level, lower_level, sep = " :: "))

# ----------------------------------------------------------------------------
# 2. THE EIGHT CLASSES AND THE LATENT STATE EACH ONE CLAIMS
# ----------------------------------------------------------------------------
# Deterministic label from three binary axes (Yhat_bin, O_local, O_rep_any):
#   1 1 1 recurrent          claims REAL       (seen everywhere)
#   1 1 0 locally unique     claims REAL       (real, locally restricted)
#   1 0 0 phantom            claims NOT-REAL   (model error, seen nowhere)
#   1 0 1 possibly missing   claims REAL       (real, undersampled locally)
#   0 0 0 possibly forbidden claims NOT-REAL   (truly absent)
#   0 0 1 locally absent     claims REAL       (real, absent under local cond.)
#   0 1 0 weak support       claims NOT-REAL   (lone sighting, likely error)
#   0 1 1 model elusive      claims REAL       (real, model missed it)
class_levels <- c("recurrent","locally unique","phantom","possibly missing",
                  "possibly forbidden","locally absent","weak support",
                  "model elusive")

# claims_real / claims_notreal : which latent state each class bets on. This is
# the bridge between the deterministic label (Framing A) and the two-state
# posterior (Framing B) — certainty is scored against the claimed state.
claims_real    <- c("recurrent","locally unique","possibly missing",
                    "locally absent","model elusive")
claims_notreal <- c("phantom","possibly forbidden","weak support")

# consistent colours (edit to match Figure 1 palette)
class_colours <- c(
  "recurrent"          = "#1b7837",
  "locally unique"     = "#7fbf7b",
  "phantom"            = "#762a83",
  "possibly missing"   = "#2166ac",
  "possibly forbidden" = "#b2182b",
  "locally absent"     = "#67a9cf",
  "weak support"       = "#f4a582",
  "model elusive"      = "#5aae61"
)

# label_from_axes : deterministic 2x2x2 lookup -> class name (Framing A)
label_from_axes <- function(yhat_bin, o_local, o_rep_any) {
  dplyr::recode(paste(yhat_bin, o_local, o_rep_any),
    "1 1 1" = "recurrent",        "1 1 0" = "locally unique",
    "1 0 0" = "phantom",          "1 0 1" = "possibly missing",
    "0 0 0" = "possibly forbidden","0 0 1" = "locally absent",
    "0 1 0" = "weak support",     "0 1 1" = "model elusive")
}

# ----------------------------------------------------------------------------
# 3. THE TWO-STATE SEQUENTIAL UPDATE (the core inference)
# ----------------------------------------------------------------------------
# prior_real : Yhat enters as the prior probability the link is real.
#              P(real) = Yhat, P(not-real) = 1 - Yhat.
#
# For each observation o in {0,1} (o=1 seen, o=0 not seen) we update P(real)
# by the likelihood ratio between the two states:
#     lr_real    = THETA_REAL^o * (1 - THETA_REAL)^(1-o)
#     lr_notreal = EPS^o        * (1 - EPS)^(1-o)
#     P(real) <- P(real)*lr_real / ( P(real)*lr_real + P(not-real)*lr_notreal )
#
# certainty : posterior mass on the state the ASSIGNED class claims.
accumulate_link <- function(yhat, o_local, rep_obs, assigned) {
  claim_is_real <- assigned %in% claims_real
  p_real <- yhat                          # prior from the model

  cert_of <- function(pr) if (claim_is_real) pr else (1 - pr)

  ev  <- c(o_local, rep_obs)              # evidence stream: local, then reps
  out <- vector("list", length(ev) + 1L)
  out[[1]] <- tibble(step = 0L, certainty = cert_of(p_real))  # step 0 = prior

  for (i in seq_along(ev)) {
    o          <- ev[i]
    lr_real    <- THETA_REAL^o * (1 - THETA_REAL)^(1 - o)
    lr_notreal <- EPS^o        * (1 - EPS)^(1 - o)
    num        <- p_real * lr_real
    p_real     <- num / (num + (1 - p_real) * lr_notreal)
    out[[i + 1L]] <- tibble(step = i, certainty = cert_of(p_real))
  }
  bind_rows(out)
}

# ----------------------------------------------------------------------------
# 4. ASSEMBLE PER-FOCAL-LINK EVIDENCE
# ----------------------------------------------------------------------------
# pair_site_gt : for each pair, the observed ground_truth at each of its sites.
#                Used to pull replicate observations (the OTHER sites) per focal.
pair_site_gt <- d %>%
  distinct(pair, focal_site, ground_truth) %>%
  group_by(pair) %>%
  summarise(sites = list(focal_site), gts = list(ground_truth), .groups = "drop")

# get_rep_obs : replicate evidence for a focal link = ground_truth at the OTHER
#               co-occurring sites of the same pair (leave-one-out).
get_rep_obs <- function(pr, fs) {
  row <- pair_site_gt[pair_site_gt$pair == pr, ]
  if (nrow(row) == 0) return(integer(0))
  sites <- row$sites[[1]]; gts <- row$gts[[1]]
  as.integer(gts[sites != fs])
}

# focal : one row per focal link, with the axes, the assigned class, and the
#         ordered replicate observations to fold in.
focal <- d %>%
  transmute(
    pair, focal_site,
    yhat     = probability,   # continuous model probability -> prior on real
    yhat_bin = prediction,    # 0/1 model call -> used for the class label
    o_local  = ground_truth   # observed at the focal site   -> O_local
  ) %>%
  mutate(
    rep_obs   = map2(pair, focal_site, get_rep_obs),
    o_rep_any = as.integer(map_int(rep_obs, ~ as.integer(any(.x == 1))) > 0),
    assigned  = label_from_axes(yhat_bin, o_local, o_rep_any),
    link_uid  = paste(pair, focal_site, sep = " @ ")
  )

# ----------------------------------------------------------------------------
# 5. RUN THE ACCUMULATION FOR EVERY LINK
# ----------------------------------------------------------------------------
# curves : link_uid x step x certainty, tagged by assigned class.
curves <- focal %>%
  mutate(traj = pmap(list(yhat, o_local, rep_obs, assigned), accumulate_link)) %>%
  select(link_uid, assigned, traj) %>%
  unnest(traj) %>%
  mutate(assigned = factor(assigned, levels = class_levels))

# ----------------------------------------------------------------------------
# 6. PLOT: per-link spaghetti, faceted by class, black median overlay
# ----------------------------------------------------------------------------
# median_curves : the class-level median trajectory (the headline shape).
median_curves <- curves %>%
  group_by(assigned, step) %>%
  summarise(certainty = median(certainty), n_links = n(), .groups = "drop")

p <- ggplot(curves, aes(step, certainty, group = link_uid, colour = assigned)) +
  geom_line(alpha = 0.06, linewidth = 0.3) +                       # one link
  geom_line(data = median_curves, aes(step, certainty, group = 1),
            colour = "black", linewidth = 0.9, inherit.aes = FALSE) +  # median
  facet_wrap(~ assigned, ncol = 4) +
  scale_colour_manual(values = class_colours, guide = "none") +
  scale_x_continuous(breaks = 0:4,
                     labels = c("prior", "local", "+1", "+2", "+3")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Evidence folded in (sequential update)",
    y = "Posterior certainty in the class assignment",
    title = "Accumulation of certainty in the taxonomic class assignment",
    subtitle = sprintf(
      "Cabrera direct observations. Illustrative rates: theta_real = %.2f, eps = %.2f. Black = median.",
      THETA_REAL, EPS)
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"))

ggsave("posterior_accumulation.pdf", p, width = 10, height = 5.5)
ggsave("posterior_accumulation.png", p, width = 10, height = 5.5, dpi = 150)

# console summary
cat("\nFocal links per class:\n")
print(focal %>% count(assigned) %>% arrange(desc(n)))
cat("\nSaved: posterior_accumulation.pdf / .png\n")
