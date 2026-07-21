# Posterior certainty in link-class assignment vs. supporting evidence
# Case study: Serra-Marin et al. (2025), Cabrera plant-pollinator direct
#             observations, with leave-one-out link predictions.
#
# The posterior over the eight taxonomic classes is
#     P(C | E)  =  P(E | C) P(C) / sum_c P(E | c) P(c)
# where E = (Yhat, O_local, O_rep).
#
# NOTE ON THE LIKELIHOOD: P(E | C) is held CONSTANT in this version (a
# placeholder until the likelihood is developed). A constant likelihood cancels
# in Bayes' rule, so here P(C | E) reduces to the prior P(C), which is shaped by
# Yhat. Consequently the posterior does not respond to supporting evidence and
# trendlines are expected to be flat. Replace `likelihood_constant()` with a
# real P(E | C) to make the x-axis informative.

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)


# ---- Parameters ----
# LIKELIHOOD_VALUE : the constant value returned for P(E | C), identical for
#                    every class. Its magnitude is irrelevant (it cancels);
#                    it exists only as the placeholder to be replaced later.
LIKELIHOOD_VALUE <- 1

# X_MAX : upper limit of the supporting-cases axis. Set to 5 as requested.
#         NOTE: under leave-one-out with 6 sites, a pair co-occurs in at most
#         4 sites, so a focal link has at most 3 other sites. Observed support
#         therefore spans 0-3 and the axis is empty beyond that.
X_MAX <- 5


# ---- Data import ----
# CSV columns (one row = one candidate species pair at one focal site):
#   higher_level : pollinator species
#   lower_level  : plant species
#   ground_truth : interaction OBSERVED at this site? (1/0)
#   probability  : model's predicted probability            -> Yhat (continuous)
#   prediction   : probability binarised to 0/1             -> Yhat_bin
#   focal_site   : which of the 6 sites is the held-out focal network
#   method       : 'obs' (direct observation) or 'rpi' (camera); we use 'obs'
df <- read.csv("results/predictions/serra_martin_loo_prediction_results.csv", stringsAsFactors = FALSE)
d  <- df %>%
  filter(method == "obs") %>%
  mutate(pair = paste(higher_level, lower_level, sep = " :: "))
# pair : unique species-pair ID, used to gather a pair's records across all
#        sites where both species co-occur (i.e. where it is a candidate link).


# ---- Evidence assembly ----
## ---- Replicate observations per pair ----
# pair_obs : for each pair, the observed ground_truth at each of its sites.
#            Used to extract the OTHER sites' evidence for each focal link.
pair_obs <- d %>%
  distinct(pair, focal_site, ground_truth) %>%
  group_by(pair) %>%
  summarise(sites = list(focal_site), obs = list(ground_truth), .groups = "drop")

# get_other_obs : replicate evidence for a focal link = the ground_truth values
#                 at the OTHER co-occurring sites of the same pair (leave-one-out).
get_other_obs <- function(pr, fs) {
  row <- pair_obs[pair_obs$pair == pr, ]
  if (nrow(row) == 0) return(integer(0))
  as.integer(row$obs[[1]][ row$sites[[1]] != fs ])
}

## ---- Per-link evidence and supporting cases ----
# yhat        : continuous model probability at the focal site (informs the prior)
# yhat_bin    : binarised model call at the focal site (defines the class label)
# o_local     : interaction observed at the focal site (1/0)
# other_obs   : vector of observations at the other co-occurring sites
# o_rep_any   : 1 if observed in at least one other site (the O_rep axis)
# support     : SUPPORTING CASES (x-axis). Number of other sites whose
#               observation AGREES with the focal prediction: if the model
#               predicted 1, a site supports by observing the interaction; if it
#               predicted 0, a site supports by not observing it.
links <- d %>%
  transmute(pair, focal_site,
            yhat     = probability,
            yhat_bin = prediction,
            o_local  = ground_truth) %>%
  mutate(
    other_obs = map2(pair, focal_site, get_other_obs),
    o_rep_any = map_int(other_obs, ~ as.integer(any(.x == 1))),
    support   = map2_int(other_obs, yhat_bin, ~ sum(.x == .y)),
    link_uid  = paste(pair, focal_site, sep = " @ ")
  )


# ---- Link classification according to taxonomy ----
# Deterministic map from the three binary axes to one of eight classes
# (Yhat_bin, O_local, O_rep_any):
#   1 1 1 recurrent           1 1 0 locally unique
#   1 0 0 phantom             1 0 1 possibly missing
#   0 0 0 possibly forbidden  0 0 1 locally absent
#   0 1 0 weak support        0 1 1 model elusive
class_levels <- c("recurrent", "locally unique", "phantom", "possibly missing",
                  "possibly forbidden", "locally absent", "weak support",
                  "model elusive")

classify <- function(yhat_bin, o_local, o_rep_any) {
  dplyr::recode(paste(yhat_bin, o_local, o_rep_any),
                "1 1 1" = "recurrent",         "1 1 0" = "locally unique",
                "1 0 0" = "phantom",           "1 0 1" = "possibly missing",
                "0 0 0" = "possibly forbidden","0 0 1" = "locally absent",
                "0 1 0" = "weak support",      "0 1 1" = "model elusive")
}

links <- links %>%
  mutate(class = factor(classify(yhat_bin, o_local, o_rep_any),
                        levels = class_levels))


# ---- Bayesian posterior over classes ----
## ---- Prior informed by the model prediction ----
# Yhat shapes which classes are plausible before consulting the field: the four
# predicted-present classes share Yhat's mass, the four predicted-absent classes
# share (1 - Yhat). Mass is spread evenly within each half (illustrative).
predicted_present <- c("recurrent", "locally unique", "phantom", "possibly missing")
predicted_absent  <- c("possibly forbidden", "locally absent", "weak support",
                       "model elusive")

make_prior <- function(yhat) {
  p <- setNames(numeric(length(class_levels)), class_levels)
  p[predicted_present] <- yhat       / length(predicted_present)
  p[predicted_absent]  <- (1 - yhat) / length(predicted_absent)
  p / sum(p)
}

## ---- Likelihood (constant placeholder) ----
# Returns P(E | C) for all eight classes. Currently constant, i.e. the evidence
# is treated as equally probable under every class. TO BE REPLACED.
likelihood_constant <- function(evidence = NULL) {
  setNames(rep(LIKELIHOOD_VALUE, length(class_levels)), class_levels)
}

## ---- Posterior of the assigned class ----
# posterior_assigned : P(C = assigned class | E), i.e. the certainty that the
#                      deterministic label given to this link is the correct one.
posterior_assigned <- function(yhat, assigned) {
  prior <- make_prior(yhat)
  lik   <- likelihood_constant()
  post  <- lik * prior            # numerator of Bayes' rule
  post  <- post / sum(post)       # normalise over the eight classes
  post[[as.character(assigned)]]
}

links <- links %>%
  mutate(posterior = map2_dbl(yhat, class, posterior_assigned))


# ---- Plotting posterior against supporting evidence ----
## ---- Palette ----
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

## ---- Scatterplot with trendlines and confidence intervals ----
# x : supporting cases (jittered slightly, as the value is a small integer)
# y : posterior certainty of the assigned class
# one colour per class, each with a linear trend and 95% CI ribbon
p <- ggplot(links, aes(x = support, y = posterior,
                       colour = class, fill = class)) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.35, size = 1.1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.15) +
  scale_colour_manual(values = class_colours, name = "Link class") +
  scale_fill_manual(values = class_colours, name = "Link class") +
  scale_x_continuous(breaks = 0:X_MAX, limits = c(-0.3, X_MAX + 0.3)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Supporting cases (sites agreeing with the focal prediction)",
    y = "Posterior certainty of the class assignment",
    title = "Certainty in link-class assignment against supporting evidence",
    subtitle = "Cabrera direct observations. Likelihood constant (placeholder): posterior reduces to the Yhat-informed prior."
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")

ggsave("posterior_vs_support.pdf", p, width = 9, height = 5.5)
ggsave("posterior_vs_support.png", p, width = 9, height = 5.5, dpi = 150)


# ---- Console summary ----
cat("\nLinks per class:\n");        print(count(links, class))
cat("\nSupporting cases observed:\n"); print(count(links, support))
cat("\nSaved: posterior_vs_support.pdf / .png\n")