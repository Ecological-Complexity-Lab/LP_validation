# ---- Bayesian reading of of predicted links taxonomy: the Cabrera worked example ----
#
# Every candidate link gets eight posterior probabilities, one per category,
# instead of a single hard label.
#
# There is one posterior function. Every SI section is a setting of its
# arguments, not a separate implementation, so the equation printed in the
# manuscript and the code producing the results are the same object.
#
# Inputs : serra_marin_loo_prediction_results.csv
# Outputs: posteriors_bayesian_example.csv (four runs x eight categories)
#          map_bayesian_richest_site.pdf

library(tidyverse)

DATA_OBS <- "results/predictions/serra_marin_loo_prediction_results.csv"
OUT_DIR  <- "results/bayesian_empirical"

# The eight categories and their error-free signatures (SI Table S2).
#   zY = the model would predict it
#   zl = it is realised at the focal site
#   zr = it is realisable somewhere in the replicates
cats <- tibble(
  category = c("recurrent", "locally unique", "possibly missing", "phantom",
               "model-elusive", "weakly-supported", "locally absent",
               "possibly forbidden"),
  zY = c(1, 1, 1, 1, 0, 0, 0, 0),
  zl = c(1, 1, 0, 0, 1, 1, 0, 0),
  zr = c(1, 0, 1, 0, 1, 0, 1, 0)
)

# Categories in which the link is neither realised locally nor realisable in
# the replicates. Feasibility confidence is one minus their posterior mass.
NOT_FEASIBLE <- c("phantom", "possibly forbidden")


# ---- 1. The evidence ----
# Each link needs four things: Y, O_l, n and R.
#
# On R. Every other site is a replicate, so R is the same for every link: the
# number of sites minus the focal one. A link can go unrecorded at a replicate
# for two reasons, and the SI does not separate them: the partners were not
# there, or they were there and the interaction was not realised or not
# detected. Both sit inside the per-replicate detection rate
# p1 = rho * (1 - eps_l) (SI Section S4), where rho is the whole ecological
# step, co-occurrence included (SI Table S7: "Co-occurrence, phenology -> rho").
# So a site where the partners were never recorded together still counts, as a
# replicate that did not record the link.
#
# There is also an ecological reason not to drop those sites. A species enters a
# site's data only if it was recorded interacting there, and co-occurring species
# are not always detected. A site without both partners in its data may well
# have held both, so excluding it would treat an undetected species as an absent
# one.
#
# n counts the replicates that recorded the link. A pair appears in a site's
# candidate grid only where both partners were recorded, so every detection
# comes from a row that exists, and a site with no row contributes n = 0.

raw <- read_csv(DATA_OBS, show_col_types = FALSE)

obs <- raw %>%
  filter(method == "obs") %>%
  mutate(pair = paste(higher_level, lower_level, sep = "||"))

N_REPLICATES <- n_distinct(obs$focal_site) - 1     # every site but the focal one

# camera records, used twice: to estimate the miss rate, and to mark the map
cam <- raw %>%
  filter(method == "rpi") %>%
  transmute(higher_level, lower_level, focal_site, camera = ground_truth)

ev <- obs %>%
  group_by(pair) %>%
  mutate(
    R = N_REPLICATES,                      # every other site
    n = sum(ground_truth) - ground_truth   # replicates that recorded the link
  ) %>%
  ungroup() %>%
  transmute(
    focal_site,
    pollinator = higher_level,
    plant      = lower_level,
    Y   = prediction,                      # hard model call
    O_l = ground_truth,                    # seen here?
    n, R,
    score = raw_score                      # raw softImpute value, see step 3
  ) %>%
  left_join(cam, by = c("pollinator" = "higher_level",
                        "plant"      = "lower_level",
                        "focal_site" = "focal_site"))

# the deterministic label, for side-by-side comparison later
ev <- ev %>%
  mutate(O_r = as.integer(n > 0)) %>%
  left_join(cats, by = character()) %>%
  filter(Y == zY, O_l == zl, O_r == zr) %>%
  select(-zY, -zl, -zr) %>%
  rename(det_category = category)


# ---- 2. The rates ----

# --- eps_Y: model error -----------------------------------------------------
# Cross-validation returns delta, the error of Y against the OBSERVATIONS, not
# eps_Y, its error against the truth-fitted prediction. No held-out set can
# give eps_Y, so delta is the stand-in, as the SI prescribes. It understates
# the model error, leaving the posterior somewhat overconfident. The SI wants a
# blocked delta; these predictions are leave-one-out by site, so this one is. 
# Here we use 1 - balanced accuracy as delta.
#
# Symmetric here, directional elsewhere, following the SI: it makes the two
# OBSERVATION axes directional (eps_l vs f, p1 vs p0) and carries the model as
# one rate. So this returns 1 - balanced accuracy and posterior() expands the
# scalar to both directions. That also keeps kappa as the SI prints it, since
# its (1-eps_Y) factor equals P(z_Y=1 | Y=1) only under symmetry. The two
# directions are far apart here (false negative 0.43, false positive 0.13);
# passing eps_Y = c(0.43, 0.13) instead raises kappa to 0.60.

estimate_eps_Y <- function(d) {
  tp <- sum(d$ground_truth == 1 & d$prediction == 1)
  fp <- sum(d$ground_truth == 0 & d$prediction == 1)
  fn <- sum(d$ground_truth == 1 & d$prediction == 0)
  tn <- sum(d$ground_truth == 0 & d$prediction == 0)
  1 - (tp / (tp + fn) + tn / (tn + fp)) / 2
}

EPS_Y <- estimate_eps_Y(obs)

# --- eps_l: local miss rate -------------------------------------------------
# From the camera contrast. Both methods watch the same site at the same time,
# so realisation is identical for both and the difference between them is
# detection alone. This is what breaks the rho / eps_l confounding: replicate
# counts alone can never separate "rarely happens" from "often missed".
#
# Conditioning on the camera-recorded links makes this DIRECTIONAL, which is
# what eps_l = P(not observed | present) needs. It is the Lincoln-Petersen
# two-sample estimator, 1 - |A n B| / |B|, and camera misses cancel from that
# ratio, so it is exactly unbiased when the two methods miss independently. A
# symmetric index is not: 1 - Jaccard also charges direct observation for the
# 124 links the CAMERA missed, giving 0.55 rather than 0.34 and dropping kappa
# from 0.533 to 0.458. Jaccard measures method overlap (J = 0.45, the two
# agreeing on under half the links they jointly recover) and belongs with the
# methodological-independence assumption, not in the likelihood.
#
# Still a lower bound, but not because the cameras miss links, which cancel.
# Detection is CORRELATED: frequent, abundant interactions are the ones both
# methods find, so the camera-recorded set is enriched for easy links and the
# miss rate measured on it understates the rate on hard ones.

estimate_eps_l <- function(obs, cam) {
  both <- obs %>% inner_join(cam, by = c("higher_level", "lower_level",
                                         "focal_site"))
  1 - sum(both$camera == 1 & both$ground_truth == 1) / sum(both$camera == 1)
}

EPS_L <- estimate_eps_l(obs, cam)

# --- f: false detection -----------------------------------------------------
# Nothing in the data identifies this. Non-co-occurring pairs, the usual
# negative control, are excluded from the candidate set by construction.
# Fixed a priori and swept in sensitivity (Section 6 below).

F_POS <- 0.05

# --- p1 and nu: replicate detection -----------------------------------------
# p1 = rho * (1 - eps_l), SI Section S4. Only the product ever enters the
# replicate factor, so we estimate it directly and never assume rho (we cannot
# separate rho from eps_l).
#
# Every site is a trial for every pair, as in Section 1, so p1 is the chance
# that a given site records the link, whatever stopped it: the partners absent,
# the interaction not realised, or the interaction missed. This is the SI's own
# estimate of p1, "the proportion of networks in which the interaction is
# recorded" (Section S9), taken over all networks rather than only those where
# the partners were recorded together. Fitting over the co-occurring sites alone
# would leave out the absent-partner step and overstate p1, and every empty
# replicate would then count as a far stronger sign of absence than it is. It
# would also rest on the recorded co-occurrence being complete, which it is not:
# co-occurring species are not always detected (Section 1).
#
# Fit to pairs recorded at least once, which is our working definition of
# "realisable in the replicates". Conditioning on n >= 1 truncates the
# distribution, so the naive proportion is biased upward. The zero-truncated
# likelihood corrects for it. The beta-binomial version also returns nu, the
# concentration that SI Section S6 needs.
#
# The homogeneous alternative is fitted alongside it. Dropping heterogeneity
# means REMOVING nu, not estimating it as infinite: the binomial is the
# beta-binomial at nu -> Inf, and under that constraint p1 is the only free
# parameter left. It has to be refitted, because the mean of a beta-binomial is
# not the binomial MLE and passing P1 with nu = Inf would be neither model.
#
# The two are nested, so which one the data prefer is testable rather than a
# matter of taste, and the likelihood ratio is reported with the rates below.

estimate_p1 <- function(obs) {
  pairs <- obs %>%
    group_by(pair) %>%
    summarise(n = sum(ground_truth), .groups = "drop") %>%
    mutate(R = N_REPLICATES + 1) %>%                  # every site is a trial
    filter(n >= 1)

  # log P(n | R) for a beta-binomial, without the binomial coefficient (it
  # cancels everywhere), minus log P(n >= 1) to correct the truncation
  nll <- function(par) {
    m <- plogis(par[1]); nu <- exp(par[2])
    a <- m * nu; b <- (1 - m) * nu
    lp <- function(k, N) lbeta(k + a, N - k + b) - lbeta(a, b)
    -sum(lp(pairs$n, pairs$R) - log1p(-exp(lp(0, pairs$R))))
  }

  # the same expression with nu removed: one detection rate shared by all links
  nll_hom <- function(par) {
    p  <- plogis(par)
    lp <- function(k, N) k * log(p) + (N - k) * log1p(-p)
    -sum(lp(pairs$n, pairs$R) - log1p(-exp(lp(0, pairs$R))))
  }

  fit <- optim(c(0, log(5)), nll)
  hom <- optimize(nll_hom, c(-8, 8))
  list(p1 = plogis(fit$par[1]), nu = exp(fit$par[2]), ll = -fit$value,
       p1_hom = plogis(hom$minimum),                  ll_hom = -hom$objective,
       n_pairs = nrow(pairs),
       n_all   = n_distinct(obs$pair),
       set     = "species pairs recorded at 1+ site, every site a trial")
}

rep_fit <- estimate_p1(obs)
P1  <- rep_fit$p1
NU  <- rep_fit$nu
P1_HOM <- rep_fit$p1_hom    # the same rate with link heterogeneity removed
P0  <- F_POS
RHO <- P1 / (1 - EPS_L)     # derived and reported, never an input

cat("\nRATES\n")
cat(sprintf("  eps_Y = %.3f   (1 - balanced accuracy)\n", EPS_Y))
cat(sprintf("  eps_l = %.3f   (camera contrast)\n", EPS_L))
cat(sprintf("  f     = %.3f   (fixed a priori)\n", F_POS))
cat(sprintf("  p1    = %.3f   nu = %.2f\n", P1, NU))
cat(sprintf("          estimated on %d of %d species pairs: %s\n",
            rep_fit$n_pairs, rep_fit$n_all, rep_fit$set))

# Heterogeneity is a nested hypothesis (nu -> Inf removes it), so the choice is
# testable. nu -> Inf sits on the boundary of the parameter space, which is why
# the null is the 50:50 mixture of chi2_0 and chi2_1 and the p-value is halved.
#
# Reading the printed line. The model with nu has one extra parameter, so it
# always fits at least as well; the question is whether the gain is worth it.
#   LR       twice the gain in log-likelihood from adding nu; 1 df = one extra
#            parameter
#   p        chance of a gain this large if all links truly shared one rate
#   dAIC     AIC with nu minus AIC without. AIC charges 2 per parameter, so
#            negative means nu still wins after paying for itself; a gap of
#            2 to 4 is weak support
LR_NU <- 2 * (rep_fit$ll - rep_fit$ll_hom)
cat(sprintf("  p1    = %.3f   with heterogeneity removed (nu absent, refitted)\n",
            P1_HOM))
cat(sprintf("          keeping nu: LR = %.2f on 1 df, boundary p = %.3f, dAIC = %.1f\n",
            LR_NU, 0.5 * pchisq(LR_NU, 1, lower.tail = FALSE), 2 - LR_NU))

# ---- 3. Calibrating the model score ----
# The `probability` column is not a probability. The prediction pipeline clips
# negative softImpute values to zero, then squashes them with a sigmoid, and a
# sigmoid maps [0, Inf) onto [0.5, 1). So the score cannot fall below 0.5 by
# construction, and a fifth of the links sit at exactly 0.5.
#
# `score` above recovers the raw value by inverting the sigmoid. Platt scaling
# is a logistic regression on that raw value, fitted site-blocked so that no
# site calibrates itself.
#
# Platt is strictly monotone, so it cannot reorder anything and the hard
# predictions must survive intact. The checklist below proves that rather than
# asserting it: the two scales give the SAME classification, and only the
# number on the threshold changes, from 0.70 on the sigmoid scale to
# Q_THRESHOLD on the calibrated one.
# Platt is monotone within a fold. Site-blocking fits six separate models, so links from different sites
# are NOT on a common scale and global monotonicity is not expected. What has
# to hold is that a single threshold still reproduces the hard predictions.

platt_calibrate <- function(ev) {
  d <- tibble(focal_site = ev$focal_site, score = ev$score, y = ev$O_l)
  q <- numeric(nrow(d))
  for (s in unique(d$focal_site)) {
    hold <- d$focal_site == s
    fit  <- glm(y ~ score, family = binomial, data = d[!hold, ])
    q[hold] <- predict(fit, newdata = d[hold, ], type = "response")
  }
  q
}

ev$q <- platt_calibrate(ev)

# the calibrated value sitting exactly at the old decision boundary
Q_THRESHOLD <- mean(c(max(ev$q[ev$Y == 0]), min(ev$q[ev$Y == 1])))

cat("\nCONSISTENCY CHECK: does Platt agree with the sigmoid classification?\n")

# Monotonicity is guaranteed within a fold, not across them, so test it there.
# A pooled Spearman below 1 is the six calibrations differing, not a failure.
fold_rho <- ev %>%
  group_by(focal_site) %>%
  summarise(rho = cor(score, q, method = "spearman"), .groups = "drop")

rank_ok <- min(fold_rho$rho)
gap_ok  <- min(ev$q[ev$Y == 1]) > max(ev$q[ev$Y == 0])
n_flip  <- sum((ev$q > Q_THRESHOLD) != (ev$Y == 1))

cat(sprintf("  [%s] Platt is monotone within every site fold (min Spearman %.4f)\n",
            ifelse(rank_ok > 0.9999, "OK", "!!"), rank_ok))
cat(sprintf("  [%s] predicted and unpredicted links do not overlap on the calibrated scale\n",
            ifelse(gap_ok, "OK", "!!")))
cat(sprintf("  [%s] thresholding q at %.4f reproduces `prediction` (%d links differ)\n",
            ifelse(n_flip == 0, "OK", "!!"), Q_THRESHOLD, n_flip))

# How far apart are the six fold-specific calibrations? Small spread means the
# score behaves the same way in rocky, dune and pine sites, so no site needs
# special handling. This is information, not a pass/fail test, and the pooled
# Spearman printed with it falls below 1 because this is six models rather than
# one, which is expected.
q_global <- predict(glm(O_l ~ score, family = binomial, data = ev),
                    type = "response")

cat(sprintf("  info calibration varies across sites by at most %.3f in q\n",
            max(abs(ev$q - q_global))))
cat(sprintf("       pooled Spearman %.4f\n",
            cor(ev$score, ev$q, method = "spearman")))

# Links the hard rule flags yet the calibration places below 0.5. This is where
# the threshold discards information: a precision-weighted rule flags them, a
# calibrated reading puts them just on the absent side. Run B treats them
# differently from Run A, which is the point of SI Fig. S9.
cat(sprintf("  note %d links are predicted (Y = 1) yet calibrate below 0.5.\n",
            sum(ev$Y == 1 & ev$q < 0.5)))

# ---- 4. The posterior master equation ----
# Calculated from the master equation the unifies all measurable elements described in the SI. 
# Nothing is hardcoded: every SI section is a setting of the arguments, not a separate branch.
#
#   Section S2  eps_r supplied, count = FALSE   symmetric binary evidence
#   Section S4  eps_r = NULL,   count = FALSE   directional rates from p1 and f
#   Section S5  count = TRUE,   nu = Inf        cumulative binomial
#   Section S6  count = TRUE,   nu finite       heterogeneous p1
#   Section S8  s = 0, pi_Y = q                 model as prior
#
# The binary forms are not an alternative formula. They are the same replicate
# distribution with the count aggregated: P(n = 0) against P(n >= 1) instead of
# P(n = k). Setting count = FALSE coarsens the evidence; it does not change the
# model. Likewise nu = Inf is the limit in which all links share one p1.
#
# eps_r is the one genuine override. Section S2 treats the replicate error as a
# free-standing rate, before it is grounded in rho and eps_l, so no value of p1
# reproduces it. Leave it NULL for everything from Section S4 onwards.
#
# Arguments
#   Y, O_l, n, R      evidence; vectors of equal length
#   eps_Y             model error; scalar, or c(false_neg, false_pos)
#   eps_l, f          local miss rate and false-detection rate
#   p1bar, nu         mean and concentration of p1 across links; nu = Inf for
#                     one shared rate
#   p0                replicate false detection; defaults to f as Table S5
#                     derives, but separable for sensitivity
#   eps_r             optional free-standing replicate rate; scalar or a pair
#   count             TRUE keeps n, FALSE collapses it to O_r = 1[n >= 1]
#   pi_Y, pi_l, pi_r  prior for each truth bit; scalars or vectors. All three
#                     at 0.5 gives 1/8 per category, the uniform prior. A bit
#                     set to 0.5 contributes the same factor to all eight
#                     categories and cancels in the normalisation, which is
#                     why pi_Y = 0.5 means the model enters only through the
#                     likelihood.
#   s                 1 = model as evidence, 0 = model as prior. SI Section S8:
#                     the two are exclusive, using both counts the model twice.
#
# Returns a matrix, one row per link, eight named columns summing to one.

posterior <- function(Y, O_l, n, R,
                      eps_Y = EPS_Y, eps_l = EPS_L, f = F_POS,
                      p1bar = P1, nu = NU, p0 = NULL,
                      eps_r = NULL, count = TRUE,
                      pi_Y = 0.5, pi_l = 0.5, pi_r = 0.5, s = 1) {

  if (is.null(p0)) p0 <- f
  if (length(eps_Y) == 1) eps_Y <- rep(eps_Y, 2)      # c(false neg, false pos)
  if (!is.null(eps_r) && length(eps_r) == 1) eps_r <- rep(eps_r, 2)

  # Beta-binomial mass without the binomial coefficient, which is identical
  # across categories and cancels. nu = Inf gives the binomial limit.
  bb <- function(k, N, m, nu) {
    if (is.infinite(nu)) return(m^k * (1 - m)^(N - k))
    a <- m * nu; b <- (1 - m) * nu
    exp(lbeta(k + a, N - k + b) - lbeta(a, b))
  }

  O_r <- as.integer(n > 0)

  L <- vapply(seq_len(nrow(cats)), function(i) {
    zY <- cats$zY[i]; zl <- cats$zl[i]; zr <- cats$zr[i]

    # --- model. The rate depends on which way the error would run.
    rY  <- if (zY == 1) eps_Y[1] else eps_Y[2]
    f_Y <- ifelse(Y == zY, 1 - rY, rY)^s

    # --- local. A realised link is missed at eps_l; an unrealised one is
    #     recorded at f.
    f_l <- if (zl == 1) ifelse(O_l == 1, 1 - eps_l, eps_l)
           else         ifelse(O_l == 1, f, 1 - f)

    # --- replicates. One distribution, read at three levels of detail.
    #     Heterogeneity applies to p1 only: p0 = f is methodological, so the
    #     non-realisable categories keep a fixed rate (Section S6).
    p_c  <- if (zr == 1) p1bar else p0
    nu_c <- if (zr == 1) nu else Inf

    f_r <- if (!is.null(eps_r)) {
      rate <- if (zr == 1) eps_r[1] else eps_r[2]
      ifelse(O_r == zr, 1 - rate, rate)
    } else if (count) {
      bb(n, R, p_c, nu_c)
    } else {
      p_none <- bb(0, R, p_c, nu_c)
      ifelse(O_r == 1, 1 - p_none, p_none)
    }

    # --- prior, Eq. S6: one probability per truth bit, multiplied.
    prior <- pi_Y^zY * (1 - pi_Y)^(1 - zY) *
             pi_l^zl * (1 - pi_l)^(1 - zl) *
             pi_r^zr * (1 - pi_r)^(1 - zr)

    f_Y * f_l * f_r * prior

  }, numeric(length(Y)))

  L <- matrix(L, nrow = length(Y), dimnames = list(NULL, cats$category))
  L / rowSums(L)
}

# feasibility confidence: realised locally, or realisable somewhere else
feasibility <- function(post) 1 - rowSums(post[, NOT_FEASIBLE, drop = FALSE])

# Maximum contextual confidence. Under asymmetric rates the SI gives kappa in
# TWO branches, because the ceiling depends on which way the local observation
# ran (SI, glossary entry for kappa):
#
#   O_l = 0   (1 - eps_Y)(1 - f)     / ((1 - f)     + eps_l)
#   O_l = 1   (1 - eps_Y)(1 - eps_l) / ((1 - eps_l) + f)
#
# They differ because a non-detection is weak evidence of absence when f is
# small, whereas a detection is strong evidence of presence. With the rates
# fitted here the two are 0.533 and 0.673, so a single number would misstate
# the ceiling for half the dataset: links recorded locally can legitimately
# rise above the O_l = 0 branch, though never above their own.
# The default is O_l = 0, which is the case the SI illustrates throughout and
# the one the accumulation figures use.
kappa <- function(eps_Y, eps_l, f, O_l = 0) {
  ifelse(O_l == 0,
         (1 - eps_Y) * (1 - f)     / ((1 - f)     + eps_l),
         (1 - eps_Y) * (1 - eps_l) / ((1 - eps_l) + f))
}


# ---- 5. Checking the function against the SI ----
# The master function is checked against the numbers printed in the SI, not against a
# second implementation. If it reproduces all of them, the manuscript and the
# code cannot drift apart.

check <- function(label, got, want, tol = 0.05) {
  cat(sprintf("  [%s] %-42s %7.1f   expected %6.1f\n",
              ifelse(abs(got - want) < tol, "OK", "!!"), label, got, want))
}

cat("\nVALIDATION AGAINST THE SI\n")

# Table S4: symmetric rates, binary evidence, uniform prior
p <- posterior(1, 0, n = 1, R = 1, eps_Y = 0.2, eps_l = 0.3, f = 0.3,
               eps_r = 0.1, count = FALSE)
check("Table S4  possibly missing", 100 * p[, "possibly missing"], 50.4)
check("Table S4  recurrent",        100 * p[, "recurrent"],        21.6)
check("Table S4  locally absent",   100 * p[, "locally absent"],   12.6)
check("Table S4  phantom",          100 * p[, "phantom"],           5.6)
check("Table S4  feasibility phi",  100 * feasibility(p),          93.0)

# Section S4: directional rates, rho = 0.15 so p1 = 0.105, R = 5
p <- posterior(1, 0, n = 1, R = 5, eps_Y = 0.2, eps_l = 0.3, f = 0.05,
               p1bar = 0.105, nu = Inf, count = FALSE)
check("Section S4  possibly missing", 100 * p[, "possibly missing"], 39.7)
check("Section S4  phantom",          100 * p[, "phantom"],          21.1)
check("Section S4  feasibility phi",  100 * feasibility(p),          73.6)

# Fig. S4a: cumulative count, detected in every replicate
for (R in c(1, 3, 5)) {
  p <- posterior(1, 0, n = R, R = R, eps_Y = 0.2, eps_l = 0.3, f = 0.05,
                 p1bar = 0.105, nu = Inf)
  check(sprintf("Fig. S4a  n = R = %d", R), 100 * p[, "possibly missing"],
        c(41.2, 54.9, 59.3)[match(R, c(1, 3, 5))])
}

# R = 0: no replicate term, so the categories pair up and split at kappa / 2
p <- posterior(1, 0, n = 0, R = 0, eps_Y = 0.2, eps_l = 0.3, f = 0.05,
               p1bar = 0.105, nu = Inf)
check("Fig. S4   R = 0 gives kappa / 2", 100 * p[, "possibly missing"],
      100 * kappa(0.2, 0.3, 0.05) / 2)
check("kappa at the SI illustrative rates", 100 * kappa(0.2, 0.3, 0.05), 60.8)

# the ceiling for a link not recorded locally: no such link can exceed it,
# however many replicates are added
cat(sprintf("\n  kappa for THIS dataset = %.3f\n", kappa(EPS_Y, EPS_L, F_POS)))


# ---- 6. The runs ----
# Two treatments of the model, crossed with two priors on the other two bits.
#
#   A  model as evidence: Y in the likelihood via eps_Y, and pi_Y = 1/2
#   B  model as prior:    pi_Y = calibrated q, and no model factor
#
# SI Section S8 forbids doing both at once, so these are separate runs.
# Comparing A with B is Fig. S9 on real data.


# --- the degree prior -------------------------------------------------------
# WHAT THIS DOES. Before looking at the evidence for a link, what should we
# believe about it? Degree answers both bits, with no free parameters.
#
#   pi_l  is it realised HERE? Under a configuration model, a pair with
#         degrees d1 and d2 in a network of L links connects at rate d1*d2/L.
#   pi_r  is it realisable ANYWHERE in the replicates? Same per-site chance,
#         R chances instead of one: 1 - (1 - pi_l)^R.
#
# There is no transfer ratio to set. The relationship between the two priors
# is derived from R, not assumed.
#
# Degrees come from the five non-focal sites, so the prior never sees the
# evidence it is about to be combined with. R is the same N_REPLICATES the
# likelihood uses (Section 1), so the prior and the evidence agree on how many
# chances a link had; co-occurrence stays in rho.

links <- obs %>% filter(ground_truth == 1) %>%
  distinct(focal_site, higher_level, lower_level)

regional <- ev %>%
  distinct(focal_site, pollinator, plant) %>%
  left_join(count(links, higher_level, name = "d_pol_all"),
            by = c("pollinator" = "higher_level")) %>%
  left_join(count(links, focal_site, higher_level, name = "d_pol_here"),
            by = c("focal_site", "pollinator" = "higher_level")) %>%
  left_join(count(links, lower_level, name = "d_pla_all"),
            by = c("plant" = "lower_level")) %>%
  left_join(count(links, focal_site, lower_level, name = "d_pla_here"),
            by = c("focal_site", "plant" = "lower_level")) %>%
  left_join(count(links, focal_site, name = "L_here"), by = "focal_site") %>%
  mutate(across(starts_with(c("d_", "L_")), ~ replace_na(.x, 0)),
         d_pol = d_pol_all - d_pol_here,
         d_pla = d_pla_all - d_pla_here,
         L_reg = nrow(links) - L_here,
         pi_l_deg = pmin(d_pol * d_pla / L_reg, 0.99),
         pi_r_deg = 1 - (1 - pi_l_deg)^N_REPLICATES) %>%
  select(focal_site, pollinator, plant, d_pol, d_pla, pi_l_deg, pi_r_deg)

ev <- ev %>% left_join(regional, by = c("focal_site", "pollinator", "plant"))

# --- the four runs ----------------------------------------------------------
# Two things can be varied independently, so we run all four combinations.
#
# AXIS 1: where the model goes.
#   A "evidence"  We use the model's hard call Y, and the likelihood pays a
#                 price eps_Y whenever Y disagrees with a category. pi_Y is
#                 then 1/2, which contributes the same factor to all eight
#                 categories and cancels in the normalisation. The model
#                 speaks once, through the likelihood.
#   B "prior"     We use the calibrated score q instead. pi_Y = q tilts the
#                 prior towards the four model-says-yes categories in
#                 proportion to how confident the model is, and the model
#                 factor is switched off in the likelihood. The model speaks
#                 once, through the prior.
#
#   These cannot be combined. SI Section S8: Y is q thresholded, so using both
#   would count the same model twice. The `s` argument enforces it.
#
#   The practical difference is the 148 links that are predicted (Y = 1) yet
#   calibrate below 0.5. Run A treats them exactly like any other predicted
#   link. Run B treats them as slightly more likely absent than present. This
#   is Fig. S9 on real data.
#
# AXIS 2: what we assume about the other two bits before seeing any evidence.
#   "uniform"  pi_l = pi_r = 1/2. Since each category is one combination of
#              the three bits, this gives (1/2)^3 = 1/8 to every category, so
#              it IS the uniform prior over the eight. Note what uniformity
#              implies: four of the eight categories have z_l = 1, so a
#              uniform prior necessarily believes a candidate link is realised
#              here with probability 1/2, against an observed rate of 0.24.
#              That is not an extra assumption, it is what uniform means here,
#              and it is a reason to run the degree prior alongside.
#   "degree"   pi_l and pi_r from the configuration model above. Two
#              generalists start higher than two specialists. Its mean pi_l is
#              0.23, essentially the observed local rate of 0.24, so on average
#              it is neither conservative nor generous. What it changes is the
#              SPREAD: specialist pairs move down and generalist pairs up,
#              where the uniform prior puts every pair at 1/2.
#
# WHAT THE COMPARISONS ANSWER.
#   A-uniform            the baseline, and the direct counterpart of the
#                        deterministic categorisation
#   A-uniform vs A-degree   how much does ecological prior information move
#                        the answer, holding the model's role fixed?
#   A-uniform vs B-uniform  how much does keeping the model's score instead of
#                        its threshold move the answer, holding the prior
#                        fixed?
#   B-degree             both together, the run using the most inputs
#
# A SECOND AXIS: THE LIKELIHOOD. The four above vary the prior and the model's
# role, all of them carrying heterogeneous p1 (Section S6). Each one also gets a
# twin with heterogeneity removed (Section S5), carrying the p1 refitted without
# it.
#
# Every prior needs its own twin, because nu and the prior do not act
# independently. nu enters the replicate factor for the four z_r = 1 categories
# only, and the posterior is normalised, so how far a change in that factor
# moves the final probabilities depends on how much prior weight sat on those
# categories to begin with. A uniform prior spreads that weight evenly; a degree
# prior concentrates it on well-connected pairs. Testing nu under the uniform
# prior alone would test it where it has the least room to matter, and the
# degree runs are the ones that carry the within-category ranking.
#
# The twins are derived from the same list rather than written out again, so a
# change to a prior cannot reach one twin and miss the other.

base_runs <- list(
  "A-uniform" = list(model_as = "evidence", pi_Y = 0.5,
                     pi_l = 0.5,         pi_r = 0.5),
  "A-degree"  = list(model_as = "evidence", pi_Y = 0.5,
                     pi_l = ev$pi_l_deg, pi_r = ev$pi_r_deg),
  "B-uniform" = list(model_as = "prior",    pi_Y = ev$q,
                     pi_l = 0.5,         pi_r = 0.5),
  "B-degree"  = list(model_as = "prior",    pi_Y = ev$q,
                     pi_l = ev$pi_l_deg, pi_r = ev$pi_r_deg)
)

# only the twins override p1bar and nu; the base runs inherit the globals
hom_runs <- setNames(
  lapply(base_runs, function(cfg) c(cfg, list(p1bar = P1_HOM, nu = Inf))),
  paste0(names(base_runs), " (no nu)"))

run_settings <- c(base_runs, hom_runs)
PRIOR_RUNS <- names(base_runs)     # the four that vary the prior
HOM_RUNS   <- names(hom_runs)      # their heterogeneity-removed twins

# The list element is called `cfg`, not `s`, because `s` is the posterior
# function's own model switch and shadowing it here would be a silent bug.
results <- imap_dfr(run_settings, function(cfg, run_name) {

  post <- posterior(
    Y     = ev$Y,
    O_l   = ev$O_l,
    n     = ev$n,
    R     = ev$R,
    eps_Y = EPS_Y,
    eps_l = EPS_L,
    f     = F_POS,
    # finite nu is heterogeneous p1 (Section S6); a run may override both to
    # remove heterogeneity, and must then supply the p1 refitted without it
    p1bar = if (is.null(cfg$p1bar)) P1 else cfg$p1bar,
    nu    = if (is.null(cfg$nu))    NU else cfg$nu,
    count = TRUE,                # keep the full count, Eq. S5
    pi_Y  = cfg$pi_Y,
    pi_l  = cfg$pi_l,
    pi_r  = cfg$pi_r,
    s     = as.integer(cfg$model_as == "evidence")
  )

  top_i <- max.col(post, ties.method = "first")

  ev %>%
    select(focal_site, pollinator, plant,
           Y, O_l, n, R, camera, det_category) %>%
    bind_cols(as_tibble(post)) %>%
    mutate(
      run          = run_name,
      model_as     = cfg$model_as,
      phi          = feasibility(post),        # feasibility confidence
      top_category = cats$category[top_i],     # category the posterior favours
      top_conf     = post[cbind(seq_len(nrow(post)), top_i)],
      det_conf     = post[cbind(seq_len(nrow(post)),
                                match(det_category, cats$category))],
      .before = 1
    )
})

write_csv(results, file.path(OUT_DIR, "posteriors_bayesian_example.csv"))

# --- headline table: hard counts against expected counts -------------------
# The deterministic version puts each link in one box; the posterior spreads
# it, so summing a column gives the expected number of links in that category.

cat("\nDETERMINISTIC COUNT vs EXPECTED COUNT (run A-uniform)\n")

one <- filter(results, run == "A-uniform") # B-degree, the most realistic model,
# gives the closest results to the deterministic categorisation

left_join(
  count(one, category = det_category, name = "deterministic"),
  one %>% select(all_of(cats$category)) %>%
    summarise(across(everything(), sum)) %>%
    pivot_longer(everything(), names_to = "category", values_to = "expected"),
  by = "category"
) %>%
  mutate(expected = round(expected, 1)) %>%
  arrange(desc(deterministic)) %>%
  print(n = 8)

# --- where the posterior disagrees with the hard label ---------------------

cat("\nWHERE THE POSTERIOR FAVOURS A DIFFERENT CATEGORY\n")

results %>%
  group_by(run) %>%
  summarise(disagreements = sum(top_category != det_category),
            mean_conf_in_hard_label = round(mean(det_conf), 3),
            .groups = "drop") %>%
  print()

# --- sensitivity to f -------------------------------------------------------
# f is the one rate with no empirical support: it is fixed a priori because
# non-co-occurring pairs, the usual negative control, are excluded from the
# candidate set by construction. So we vary it and report how much moves.
# f enters three places: the local false-positive factor, p0 on the replicate
# axis, and kappa.

sensitivity_f <- map_dfr(seq(0, 0.15, by = 0.01), function(f_try) {
  post <- posterior(ev$Y, ev$O_l, ev$n, ev$R, f = f_try, p0 = f_try)
  tibble(
    f                = f_try,
    kappa            = kappa(EPS_Y, EPS_L, f_try),
    expected_missing = sum(post[, "possibly missing"]),
    expected_phantom = sum(post[, "phantom"]),
    mean_phi         = mean(feasibility(post))
  )
})

cat("\nSENSITIVITY TO f\n")
print(sensitivity_f, n = 16)

# The number worth quoting: how far the headline moves across the whole range.
cat(sprintf("\n  Expected 'possibly missing' ranges %.1f to %.1f as f goes 0 to 0.15.\n",
            min(sensitivity_f$expected_missing),
            max(sensitivity_f$expected_missing)))

# ---- 7. The map ----
# The richest site, one category, shaded by posterior. Same species, same
# ordering and same look as map_richest_site in the deterministic script, so
# the two can be laid side by side. An asterisk marks links the cameras found.

# -- Top-of-scale colour per category, from map_colors in
#    serra_marin_link_classification_clean.R (Section 2). White is the floor
#    everywhere, so each map runs white -> the colour that category carries in
#    the deterministic figure, and the two are read with one key. The published
#    palette splits phantom and possibly forbidden by whether the other
#    sampling method corroborated them; there is no such split here, so both
#    take their no-evidence colour. --
CATEGORY_COLOUR <- c(
  "recurrent"          = "#C05030",   # dark coral
  "model-elusive"      = "#8A3050",   # deep rose
  "possibly missing"   = "#DBA040",   # amber
  "locally absent"     = "#F9C8B0",   # pale peach
  "locally unique"     = "#B0C4DE",   # light blue
  "weakly-supported"   = "#E8A8B8",   # pink
  "phantom"            = "#D5C8F0",   # light lavender  (no-evidence tint)
  "possibly forbidden" = "#5E3DA0"    # dark purple     (no-evidence tint)
)

# Richest site by SPECIES richness (plants + pollinators), which is the
# criterion the deterministic script uses. Counting links instead happens to
# pick the same site here, but the two are different questions.
richest <- obs %>%
  filter(ground_truth == 1) %>%
  group_by(focal_site) %>%
  summarise(n_species = n_distinct(lower_level) + n_distinct(higher_level),
            .groups = "drop") %>%
  slice_max(n_species, n = 1) %>%
  pull(focal_site)

# Species ordered by overall degree, matching map_rich in the deterministic
# script. Degree is the number of DISTINCT partners pooled across sites, which
# is not the same as counting (site, partner) rows: the latter sums the
# per-site degrees and gives a different order.
deg_pol <- obs %>% filter(ground_truth == 1) %>%
  group_by(higher_level) %>%
  summarise(deg = n_distinct(lower_level), .groups = "drop")
deg_pla <- obs %>% filter(ground_truth == 1) %>%
  group_by(lower_level) %>%
  summarise(deg = n_distinct(higher_level), .groups = "drop")

ital <- function(x) lapply(strsplit(x, " "),
                           function(y) bquote(italic(.(paste(y, collapse = " ")))))

# make_map(run, category)
#
# One map: the chosen site, shaded by the posterior for the chosen category
# under the chosen run. Geometry, typography and species ordering are taken
# from map_rich in the deterministic script, so any two of these panels, or one
# of these and the deterministic map, can sit side by side without rescaling.
#
#   run       one of the names in run_settings ("A-uniform", "A-degree",
#             "B-uniform", "B-degree"). A-uniform is the direct counterpart of
#             the deterministic map.
#   category  one of the eight in cats$category. Sets both the shaded quantity
#             and the top-of-scale colour.
#   site      defaults to the richest site.
#   camera_marks  which links get an asterisk for a camera record.
#             "assigned"  only links the deterministic rule put in THIS
#                         category. The asterisks then read as corroboration of
#                         that category: of the links called possibly missing,
#                         these are the ones the independent method found.
#             "all"       every camera record at the site, whatever category it
#                         was assigned. Use this to see camera coverage as a
#                         backdrop rather than as a verdict on one category.
#             "none"      no asterisks.
#   frame_assigned  outline in black the cells the DETERMINISTIC rule put in
#             this category. Colour is the posterior, spread over every cell,
#             so nothing otherwise shows where the hard rule drew its line.
#             With the frames both readings are visible: a pale framed cell is
#             a link the rule assigned but the posterior doubts, a dark
#             unframed cell one the posterior favours but the rule gave away.
#   limit     top of the colour scale. NULL stops it at the highest posterior
#             actually reached, which is the honest default: kappa caps every
#             posterior, so a 0-to-1 scale would wash the figure out and
#             overstate confidence. Pass a number to hold several maps on one
#             scale, which is what makes them comparable to each other.
#
# Returns the ggplot, with the scale top attached as attr(, "conf_max"), the
# asterisk count as attr(, "n_marked") and the framed count as attr(, "n_framed").
make_map <- function(run      = "A-uniform",
                     category = "possibly missing",
                     site     = richest,
                     results_df = results,
                     camera_marks = c("assigned", "all", "none"),
                     frame_assigned = TRUE,
                     limit    = NULL,
                     # caption = FALSE drops the footnote so the grid itself
                     # renders larger. Use it when the surrounding document
                     # already states the site, scale, frames and marks.
                     caption  = FALSE) {

  camera_marks <- match.arg(camera_marks)

  if (!category %in% cats$category)
    stop("category must be one of: ", paste(cats$category, collapse = ", "))
  if (!run %in% unique(results_df$run))
    stop("run must be one of: ", paste(unique(results_df$run), collapse = ", "))

  d <- results_df %>%
    filter(run == !!run, focal_site == !!site) %>%
    mutate(conf = .data[[category]])

  if (nrow(d) == 0) stop("no rows for run '", run, "' at site '", site, "'")

  # Levels are taken from this site's own rows, because several species tie on
  # degree and a tie is broken by the order the rows arrive in. Ordering
  # globally and then subsetting picks a different winner and shifts the axis.
  pol_levels <- d %>% distinct(pollinator) %>%
    left_join(deg_pol, by = c("pollinator" = "higher_level")) %>%
    arrange(desc(deg)) %>% pull(pollinator)
  pla_levels <- d %>% distinct(plant) %>%
    left_join(deg_pla, by = c("plant" = "lower_level")) %>%
    arrange(desc(deg)) %>% pull(plant)

  # Not reversed on the y axis: ggplot puts the first level at the bottom,
  # which is where the deterministic map puts the highest-degree plant.
  d <- d %>%
    mutate(pollinator = factor(pollinator, levels = pol_levels),
           plant      = factor(plant,      levels = pla_levels))

  conf_max <- if (is.null(limit)) max(d$conf) else limit

  # Which links carry an asterisk. camera == 1 drops NA of its own accord, so
  # links with no camera coverage are simply unmarked, the same as links the
  # cameras watched and did not find. The map cannot tell those two apart; the
  # corroboration table in Section 6 is where that distinction lives.
  marked <- switch(camera_marks,
    none     = d[0, ],
    all      = filter(d, camera == 1),
    assigned = filter(d, camera == 1, det_category == category)
  )

  cap_marks <- switch(camera_marks,
    none     = "",
    all      = sprintf(". * marks a camera record (%d)", nrow(marked)),
    assigned = sprintf(". * marks a camera record among the %d links assigned to this category (%d of them)",
                       sum(d$det_category == category), nrow(marked))
  )

  # cells the deterministic rule assigned to this category
  framed <- if (frame_assigned) filter(d, det_category == category) else d[0, ]
  cap_frame <- if (nrow(framed) > 0)
    sprintf(". Black outlines are the %d links the deterministic rule assigned here",
            nrow(framed)) else ""

  p <- ggplot(d, aes(pollinator, plant, fill = conf)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = unname(CATEGORY_COLOUR[category]),
                        limits = c(0, conf_max),
                        name = sprintf("P(%s)", category)) +
    scale_x_discrete(labels = ital) +
    scale_y_discrete(labels = ital) +
    labs(x = "Pollinator", y = "Plant",
         caption = if (caption)
           sprintf("%s, run %s. Scale stops at %.2f%s%s.",
                   site, run, conf_max, cap_frame, cap_marks) else NULL) +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(size = 16, angle = 90, vjust = 0.5),
      axis.text.y  = element_text(size = 16),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 14),
      panel.grid   = element_blank()
    )

  # Drawn before the asterisks so a mark is never hidden by a frame edge.
  # inherit.aes = FALSE keeps this layer off the fill scale, so the outline is
  # added without a second entry appearing in the legend.
  if (nrow(framed) > 0) {
    p <- p + geom_tile(data = framed, aes(pollinator, plant),
                       inherit.aes = FALSE, fill = NA,
                       colour = "black", linewidth = 0.45)
  }

  if (nrow(marked) > 0) {
    p <- p + geom_text(data = marked, aes(label = "*"),
                       colour = "grey20", size = 5, vjust = 0.72)
  }

  attr(p, "conf_max") <- conf_max
  attr(p, "n_marked") <- nrow(marked)
  attr(p, "n_framed") <- nrow(framed)
  p
}

# --- what the asterisks mark, by category -----------------------------------
# The two marking modes side by side. camera_marks = "all" puts the same
# `marks_all` asterisks on every map, whatever category is shaded;
# camera_marks = "assigned" splits that same total across the categories, so
# the `confirmed` column sums to it. That is the corroboration result: of the
# links the deterministic rule put in each category, how many did the
# independent method find?
#
# `covered` is the column that keeps the rates honest. A link with no camera
# coverage is unmarked for the same reason as a link the cameras watched and
# did not find, and the map cannot tell them apart. rate_covered divides by
# the links the cameras actually watched; rate_assigned divides by all of them
# and therefore understates corroboration wherever coverage is thin.

camera_mark_summary <- function(run = "A-uniform", site = richest,
                                results_df = results) {

  d <- results_df %>% filter(run == !!run, focal_site == !!site)
  if (nrow(d) == 0) stop("no rows for run '", run, "' at site '", site, "'")

  d %>%
    group_by(category = det_category) %>%
    summarise(
      assigned  = n(),
      covered   = sum(!is.na(camera)),
      confirmed = sum(camera == 1, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    right_join(tibble(category = cats$category), by = "category") %>%
    mutate(across(c(assigned, covered, confirmed), ~ replace_na(.x, 0)),
           rate_assigned = confirmed / pmax(assigned, 1),
           rate_covered  = ifelse(covered > 0, confirmed / covered, NA_real_)) %>%
    arrange(desc(rate_covered), desc(confirmed))
}

MARK_RUN <- "A-uniform"
mark_tbl <- camera_mark_summary(MARK_RUN)

# The total below is what camera_marks = "all" puts on every map. The
# "assigned" mode splits that same total across the categories, so `confirmed`
# is what each map marks and the column sums back to it.
cat(sprintf("\nCAMERA CORROBORATION AT %s (run %s): %d records\n",
            richest, MARK_RUN, sum(mark_tbl$confirmed)))

mark_tbl %>%
  mutate(rate_assigned = sprintf("%.1f%%", 100 * rate_assigned),
         rate_covered  = ifelse(is.na(rate_covered), "-",
                                sprintf("%.1f%%", 100 * rate_covered))) %>%
  print(n = 8, width = Inf)

# --- produce the maps -------------------------------------------------------
# Three categories x four runs x three formats = 36 files. Written as a grid
# rather than as repeated blocks so a change to size or marking applies to all
# of them at once, and so no combination can be silently missed or misnamed.
#
# Each map is also left in the session under the name it is saved as, e.g.
# map_bayes_forbidden_A_uniform, for interactive use.

MAP_CATEGORIES <- c("possibly forbidden", "possibly missing", "phantom")
# single maps for the four prior runs. The heterogeneity twins are shown as
# their own four-panel figures instead, which is where they are read against
# the four-panel figure they mirror.
MAP_RUNS       <- PRIOR_RUNS
MAP_FORMATS    <- c("pdf", "svg", "png")
MAP_MARKS      <- "assigned"     # asterisks only on this category's own links
MAP_W          <- 14
MAP_H          <- 7

# short, file-safe stems: category words, and the run name with anything not
# valid in a filename or an object name collapsed to "_"
map_stem <- c("possibly forbidden" = "forbidden",
              "possibly missing"   = "missing",
              "phantom"            = "phantom")

run_stem <- function(x) gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", x))

map_grid <- expand_grid(category = MAP_CATEGORIES, run = MAP_RUNS) %>%
  mutate(object = sprintf("map_bayes_%s_%s",
                          map_stem[category], run_stem(run)))

map_log <- vector("list", nrow(map_grid))

for (i in seq_len(nrow(map_grid))) {
  cat_i <- map_grid$category[i]
  run_i <- map_grid$run[i]
  obj_i <- map_grid$object[i]

  p <- make_map(run = run_i, category = cat_i, camera_marks = MAP_MARKS)
  assign(obj_i, p)                       # available in the session by name

  for (fmt in MAP_FORMATS) {
    ggsave(file.path(OUT_DIR, paste0(obj_i, ".", fmt)), p,
           width = MAP_W, height = MAP_H)
  }

  map_log[[i]] <- tibble(
    object   = obj_i,
    run      = run_i,
    category = cat_i,
    conf_max = attr(p, "conf_max"),
    asterisks = attr(p, "n_marked")
  )
}

map_log <- bind_rows(map_log)

cat(sprintf("\n\nMAPS WRITTEN: %d combinations x %d formats = %d files in %s\n",
            nrow(map_grid), length(MAP_FORMATS),
            nrow(map_grid) * length(MAP_FORMATS), OUT_DIR))
cat(sprintf("site: %s, camera marks: %s\n\n", richest, MAP_MARKS))
# conf_max in the log is where each colour scale tops out, so it differs from
# map to map. Pass limit = <value> to make_map() to hold several on one scale.
map_log %>%
  mutate(conf_max = round(conf_max, 3)) %>%
  arrange(category, match(run, MAP_RUNS)) %>%
  print(n = Inf, width = Inf)

# --- the four runs side by side, for one category ---------------------------
# make_run_panel("phantom") builds one figure holding all four runs, so the
# effect of the prior is read off a single image rather than by flicking
# between files.
#
# The important argument is `shared_limit`. Each map on its own stretches the
# colour ramp to its own maximum, and those maxima differ a lot between runs,
# most of all between the uniform and the degree priors.
# Four panels on four different scales would make an identical shade mean four
# different probabilities, so the default holds every panel on one scale, the
# highest reached by any of them. That also lets the four legends collapse into
# one. Set shared_limit = FALSE only to inspect a single run's internal
# structure, never to compare runs.
#
# Axis labels are drawn once per edge rather than in all four panels: the
# species are identical everywhere, so repeating 46 names four times spends
# most of the figure on text.
#
#   category      one of cats$category
#   runs          which runs, in order; defaults to all four
#   site, camera_marks, frame_assigned, ...  passed through to make_map().
#                 The black frames are identical in every panel, since the
#                 deterministic assignment does not depend on the run, so they
#                 double as a fixed reference against which the shading moves.
#   ncol          2 gives a 2 x 2 block; 1 or 4 give a strip
#   file          output path; NULL uses a name built from the category
#
# Returns the patchwork object, with the shared scale top as attr(, "limit").

library(patchwork)

make_run_panel <- function(category,
                           # the four prior/model runs. The likelihood
                           # sensitivity run gets its own single maps instead,
                           # so this figure stays a 2 x 2 comparison of priors.
                           runs         = PRIOR_RUNS,
                           site         = richest,
                           camera_marks = "assigned",
                           frame_assigned = TRUE,
                           shared_limit = TRUE,
                           # FALSE drops the shared footnote, which gives the
                           # four grids the space it was using
                           caption      = FALSE,
                           ncol         = 2,
                           width = 24, height = 13, dpi = 150,
                           file = NULL,
                           save = TRUE) {

  if (!category %in% cats$category)
    stop("category must be one of: ", paste(cats$category, collapse = ", "))

  # one scale for all panels: the highest posterior any run reaches here
  lim <- if (shared_limit) {
    max(vapply(runs, function(r)
      max(results[[category]][results$run == r & results$focal_site == site]),
      numeric(1)))
  } else NULL

  raw_maps <- lapply(runs, function(r)
    make_map(run = r, category = category, site = site,
             camera_marks = camera_marks, frame_assigned = frame_assigned,
             limit = lim))
  # read attributes BEFORE adding to the plots, since ggplot drops them
  n_marked <- attr(raw_maps[[1]], "n_marked")
  n_framed <- attr(raw_maps[[1]], "n_framed")

  panels <- lapply(seq_along(runs), function(i) {
    p <- raw_maps[[i]]
    # the shared caption below carries site and scale, so drop the per-panel one
    p <- p + labs(title = runs[i], caption = NULL) +
      theme(plot.title = element_text(size = 15, face = "bold", hjust = 0))
    # keep axis text only on the outer edges
    bottom_row <- i > (length(runs) - ncol)
    left_col   <- (i - 1) %% ncol == 0
    if (!bottom_row) p <- p + theme(axis.text.x  = element_blank(),
                                    axis.title.x = element_blank(),
                                    axis.ticks.x = element_blank())
    if (!left_col)   p <- p + theme(axis.text.y  = element_blank(),
                                    axis.title.y = element_blank(),
                                    axis.ticks.y = element_blank())
    p
  })

  cap <- paste0(
    sprintf("%s. Colour is P(%s); all panels share one scale topping at %.2f, so shades are comparable between runs.",
            site, category, if (is.null(lim)) NA_real_ else lim),
    # the frames are the same in every panel: the hard rule does not vary by run
    if (n_framed > 0)
      sprintf(" Black outlines are the %d links the deterministic rule assigned here, identical in all panels.",
              n_framed) else "",
    if (identical(camera_marks, "none")) ""
    else sprintf(" * marks a camera record (%d).", n_marked))

  fig <- wrap_plots(panels, ncol = ncol, guides = "collect") +
    plot_annotation(tag_levels = "a",
                    caption = if (caption) cap else NULL) &
    theme(plot.tag = element_text(size = 16, face = "bold"))

  if (save) {
    if (is.null(file))
      file <- file.path(OUT_DIR,
                        sprintf("panel_runs_%s.png",
                                if (category %in% names(map_stem))
                                  map_stem[[category]]
                                else gsub("[^a-z]+", "_", category)))
    ggsave(file, fig, width = width, height = height, dpi = dpi, bg = "white")
    cat(sprintf("  panel written: %s   (shared scale top %.3f)\n",
                file, if (is.null(lim)) NA_real_ else lim))
  }

  attr(fig, "limit") <- lim
  fig
}

# caption = FALSE: these panels go into note_empirical_bayesian_overview.Rmd,
# which states the site, scale, frames and marks in its own prose, so the
# footnote is redundant there and the four grids get the space instead.
#
# Two figures per category: the four prior runs, and their heterogeneity twins
# in the same panel order, so the pair can be read side by side. Each figure
# takes its own shared colour scale, so compare shapes across the pair and
# shades only within a figure.
#
# Figures are printed only in an interactive session. Under Rscript, R would
# send them to its default device and leave a stray Rplots.pdf duplicating the
# files already written.
cat("\n\nFOUR-RUN PANELS\n")
for (k in MAP_CATEGORIES) {
  panel_het <- make_run_panel(k, caption = FALSE)
  panel_hom <- make_run_panel(
    k, runs = HOM_RUNS, caption = FALSE,
    file = file.path(OUT_DIR, sprintf("panel_runs_%s_no_nu.png", map_stem[k])))
  if (interactive()) {
    print(panel_het)        # the four prior runs
    print(panel_hom)        # their heterogeneity twins
  }
}

# Any other combination is one call, for example
#   make_map("B-degree", "phantom")
#   make_map("A-uniform", "possibly forbidden")
# every camera record at the site rather than only this category's:
#   make_map("A-uniform", "phantom", camera_marks = "all")
# and a shared `limit` puts several on one scale:
#   make_map("A-uniform", "possibly missing", limit = 0.6)


# ---- 8. Category distributions across the runs ----
# Two views of the same posteriors, sharing one category colour key taken from
# the maps. Produced twice: once for the four prior/model runs, and once for
# A-uniform against A-uniform (no nu), which isolates the likelihood.
#
#   a  where the mass goes. Summing a posterior column gives the EXPECTED
#      number of links in that category, which is the honest headline number
#      once the hard label is replaced by a distribution. The tick marks the
#      deterministic count, identical in all four panels, so the gap between
#      tick and bar is exactly what each prior redistributes.
#   b  how confident the posterior gets. The spread of P(category) across every
#      link, which says whether a category's mass is a few confident links or
#      many diffuse ones.
#
# On the kappa line in panel b. kappa, the ceiling on confidence, has two
# versions depending on what the local observation said (see kappa() in
# Section 4):
#
#   O_l = 0, link not recorded locally  (1 - eps_Y)(1 - f)     / ((1 - f) + eps_l)
#   O_l = 1, link recorded locally      (1 - eps_Y)(1 - eps_l) / ((1 - eps_l) + f)
#
# The line in panel b is the first. It carries TWO conditions, and points sit to
# the right of it whenever either is relaxed, which is a result rather than a
# violation.
#
#   1. A uniform prior. B replaces the flat pi_Y with the calibrated score and
#      the degree runs add ecological priors, so the line is drawn only in the
#      uniform-prior panels; the others exceed it and should.
#   2. A link NOT recorded locally (O_l = 0). The derivation is the ceiling for
#      a non-detection, where the local factor contributes only eps_l against
#      1 - f. A positive local detection contributes 1 - eps_l against f, a far
#      sharper contrast, so those links can go higher, up to the O_l = 1
#      version of kappa above. The ones that do are recurrent or model-elusive,
#      the two categories with z_l = 1 and z_r = 1.
#
# Each row mixes both local conditions, so the line cannot be drawn per row.
# The axis label states what it bounds instead.

CAT_ORDER <- cats$category          # taxonomy order: predicted four, then not

# the deterministic counts are a property of the labelling, not of a run, so
# the same reference is drawn in every panel
det_counts <- results %>%
  filter(run == names(run_settings)[1]) %>%
  count(category = det_category, name = "deterministic") %>%
  mutate(category = factor(category, levels = rev(CAT_ORDER)))

si_base <- theme_classic(base_size = 10) +
  theme(axis.line        = element_line(colour = "grey40", linewidth = 0.3),
        axis.ticks       = element_line(colour = "grey40", linewidth = 0.3),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 9, hjust = 0),
        panel.spacing    = unit(9, "pt"),
        plot.tag         = element_text(size = 13, face = "bold"),
        legend.position  = "none")

# Built as a function so the same two panels serve both the prior comparison
# and the heterogeneity sensitivity, rather than the second being a copy that
# can drift from the first.
#   runs        which runs to show, left to right
#   kappa_runs  panels that get the kappa line: uniform-prior runs only
distribution_figure <- function(runs, stem, width = 11, height = 7,
                                kappa_runs = runs[1]) {

  post_long <- results %>%
    filter(run %in% runs) %>%
    select(run, all_of(CAT_ORDER)) %>%
    pivot_longer(all_of(CAT_ORDER), names_to = "category", values_to = "p") %>%
    mutate(category = factor(category, levels = rev(CAT_ORDER)),
           run      = factor(run, levels = runs))

  expected_tbl <- post_long %>%
    group_by(run, category) %>%
    summarise(expected = sum(p), .groups = "drop") %>%
    left_join(det_counts, by = "category")

  # kappa does not involve p1 or nu, so it bounds the homogeneous run exactly
  # as it bounds A-uniform
  kappa_line <- tibble(run = factor(kappa_runs, levels = runs),
                       k   = kappa(EPS_Y, EPS_L, F_POS))

  p_mass <- ggplot(expected_tbl, aes(expected, category, fill = category)) +
    geom_col(width = 0.72, colour = "grey35", linewidth = 0.2) +
    geom_point(aes(x = deterministic), shape = 124, size = 2.6,
               colour = "grey15") +
    facet_wrap(~ run, nrow = 1) +
    scale_fill_manual(values = CATEGORY_COLOUR) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
    labs(x = "Expected number of links   (| = deterministic count)", y = NULL) +
    si_base

  p_spread <- ggplot(post_long, aes(p, category, fill = category)) +
    geom_vline(data = kappa_line, aes(xintercept = k), linetype = "22",
               colour = "grey45", linewidth = 0.35) +
    # outliers carry the story here (the links that pass kappa), so they are
    # drawn large enough to count rather than as faint dust
    geom_boxplot(width = 0.62, colour = "grey35", linewidth = 0.25,
                 outlier.size = 0.55, outlier.colour = "grey30",
                 outlier.alpha = 0.45, outlier.stroke = 0) +
    geom_text(data = kappa_line, aes(x = k, y = 8.72, label = "kappa"),
              parse = TRUE, inherit.aes = FALSE, hjust = -0.2, size = 2.9,
              colour = "grey35") +
    facet_wrap(~ run, nrow = 1) +
    scale_fill_manual(values = CATEGORY_COLOUR) +
    # headroom so the kappa label sits inside the panel rather than clipped
    scale_y_discrete(expand = expansion(add = c(0.6, 1.0))) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = expression("P(category) per link"), y = NULL) +
    si_base +
    theme(strip.text = element_blank())    # run names already label panel a

  fig <- (p_mass / p_spread) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 13, face = "bold"))

  ggsave(file.path(OUT_DIR, paste0(stem, ".pdf")), fig,
         width = width, height = height)
  ggsave(file.path(OUT_DIR, paste0(stem, ".png")), fig,
         width = width, height = height, dpi = 300)
  invisible(list(fig = fig, expected = expected_tbl))
}

main_dist <- distribution_figure(PRIOR_RUNS, "fig_category_distributions")
if (interactive()) print(main_dist$fig)

cat("\nEXPECTED COUNTS BY RUN (deterministic count in the first column)\n")
main_dist$expected %>%
  mutate(expected = round(expected, 1)) %>%
  select(category, run, deterministic, expected) %>%
  pivot_wider(names_from = run, values_from = expected) %>%
  arrange(desc(deterministic)) %>%
  print(n = 8, width = Inf)

# --- heterogeneity sensitivity, as its own figure ---------------------------
# The same figure for the four twins, in the same panel order, so it can be
# read directly against the one above. Each pair of panels differs by the
# likelihood alone: same prior, same model role, same rates, one carrying
# heterogeneous p1 and the other the single p1 refitted without it. Anything
# that moves between a panel and its counterpart is nu, and nothing else.
hom_dist <- distribution_figure(HOM_RUNS, "fig_category_distributions_no_nu",
                                kappa_runs = HOM_RUNS[1])
if (interactive()) print(hom_dist$fig)

# Paired by prior, since the question is what nu does WITHIN a prior, not which
# prior is closer to the deterministic counts.
cat("\nHETEROGENEITY SENSITIVITY: shift in expected counts when nu is removed\n")
bind_rows(main_dist$expected %>% mutate(nu = "het"),
          hom_dist$expected  %>% mutate(nu = "hom",
                                        run = sub(" \\(no nu\\)$", "", run))) %>%
  mutate(run = as.character(run)) %>%
  select(category, run, nu, deterministic, expected) %>%
  pivot_wider(names_from = nu, values_from = expected) %>%
  mutate(shift = round(hom - het, 1)) %>%
  select(category, run, deterministic, shift) %>%
  pivot_wider(names_from = run, values_from = shift) %>%
  arrange(desc(deterministic)) %>%
  print(n = 8, width = Inf)


# ---- 9. Camera validation: does P(category) predict corroboration? ----
#
# RATIONALE
# The posterior is built from three sources: the model, the direct observation
# at the focal site, and the replicate counts. The cameras are a fourth source
# it never saw. Each category makes a claim about what is happening at the
# focal site, and each claim implies whether an independent method should or
# should not find the link. Testing those implications validates the
# CATEGORIES rather than the arithmetic that produced them.
#
# THE TEST
# One logistic regression per category, over every link a camera watched:
#
#     camera ~ P(category)                       raw
#     camera ~ P(category) + O_l                 direct observation held fixed
#     camera ~ P(category) + O_l + focal_site    and camera placement too
#
# reported as an odds ratio per one standard deviation of that category's
# posterior. The raw model is close to tautological, because the categories are
# defined partly BY O_l, so P(category) largely re-expresses the observation.
# The second model is the real test: among links whose direct observation was
# the same, does the posterior still say which ones the cameras find?
#
# VARIABLES
#   camera        outcome. Did the cameras record this pair at this site? 1/0,
#                 and NA where the pair falls outside the camera grid.
#                 What NA means: the candidate set for a method at a site is the
#                 full plant x pollinator grid of the species THAT METHOD
#                 recorded there, so a link has a camera row only when BOTH
#                 partners appear in that site's camera data. Of the 673 links
#                 without one, 424 lack the pollinator, 158 the plant, 91 both.
#                 So camera = 0 means the cameras caught both species at this
#                 site but never together, which is a real negative; it does NOT
#                 mean a camera was trained on that particular pair. Those NA
#                 links are dropped, since they can neither confirm nor refute.
#                 Note the selection this implies: a pollinator enters a site's
#                 camera set only if the cameras detected it there at all, so
#                 the graded subset is enriched for camera-detectable species
#                 and coverage ranges from 41% to 82% across sites. That is one
#                 reason focal_site is carried as a covariate below.
#   P(category)   predictor. The posterior for the category under test,
#                 standardised so the odds ratio is per SD and categories with
#                 different spreads stay comparable.
#   O_l           covariate. The direct-observation record, 1/0.
#   focal_site    covariate. Six levels. Cameras were not placed at random, so
#                 a category could look corroborated merely by sitting at
#                 well-watched sites.
#
# PARAMETERS
#   None are fitted here beyond the regression coefficients. This test needs no
#   camera detection rate and no false-positive rate, unlike a calibration
#   check: it asks only whether the posterior ORDERS links correctly, not
#   whether its absolute values are right. The posterior itself still carries
#   eps_l, which was estimated from these same cameras, so the evidence is not
#   fully independent. Holding O_l fixed is what limits that exposure, since
#   eps_l enters through the local axis.
#
# READING THE OUTPUT
#   n, confirmed   size of the DETERMINISTIC category AMONG LINKS INSIDE THE
#                  CAMERA GRID, and how many of those the cameras recorded. Not
#                  the size of the whole category: possibly forbidden has 941
#                  links in total and 452 here. Descriptive only: every
#                  regression below
#                  uses all links inside the camera grid, not just these.
#   OR_adj, p_adj  the decisive columns. Odds ratio per 1 SD of this category's
#                  posterior, with O_l held fixed, and its p-value.
#                  OR > 1  a higher posterior goes with MORE camera records
#                  OR < 1  a higher posterior goes with FEWER camera records
#                  OR = 1  the posterior says nothing the observation had not.
#                  An odds ratio is a multiplier on the odds p/(1-p), not on the
#                  probability, so 1.6 is NOT "1.6 times as likely".
#   pct_avg        the same effect on the probability scale, which is what to
#   pct_1sd        quote in prose. Every link is predicted twice from the
#   pct_diff       O_l-adjusted model, once at the average posterior (pct_avg)
#                  and once one SD above it (pct_1sd), then averaged across
#                  links; pct_diff is the gap in percentage points. Averaging
#                  is necessary because a fixed odds ratio implies DIFFERENT
#                  probability changes at different starting points, smaller
#                  among locally unobserved links, where camera records are
#                  rare, than among observed ones. Averaging over the links
#                  respects the real mix of O_l and answers "for a typical link
#                  here, how much does the chance of a camera record move?".
#                  Read these together with p_adj: a tiny p with a gap of a few
#                  points is a clean signal in a large sample, not a large
#                  effect.
#   verdict        does the significant effect run the way the category
#                  predicts? The expected direction is set below per category,
#                  from what the category claims ecologically rather than from
#                  its z_l bit. Those differ: possibly missing carries z_l = 0,
#                  yet its whole claim is that the interaction IS present and
#                  was missed, so the cameras should find it MORE often.
#   OR_site        as OR_adj but also adjusting for focal_site, since camera
#   p_site         coverage runs from 41% to 82% across sites.
#   OR_raw, p_raw  reference only, no covariates. Close to tautological, since
#                  categories are defined partly BY O_l. It can even carry the
#                  opposite sign to OR_adj, as possibly missing has done. Never
#                  read it against the verdict.

# What each category claims, and therefore what the cameras should show.
cam_expectation <- tibble::tribble(
  ~category,             ~expect,     ~rationale,
  "recurrent",           "higher",
  "Predicted and recorded both here and in replicates. Every source agrees the interaction is real, so an independent method should find it too.",
  "locally unique",      "higher",
  "Predicted and recorded here but in no replicate. It is real at this site, and the absence elsewhere says nothing about here, so the camera should find it.",
  "possibly missing",    "higher",
  "Predicted and recorded in replicates but not here. The claim is that the interaction is present and direct observation missed it, so an independent method should recover it more often than for comparable unobserved links.",
  "phantom",             "unclear",
  "Predicted but recorded nowhere. Either a model false positive, in which case the camera finds nothing, or an interaction this method cannot see, in which case it does. The two readings predict opposite outcomes, which is why the prescribed action is an independent method rather than more effort.",
  "model-elusive",       "higher",
  "Recorded here and in replicates but never predicted. The observations already agree it is real and the failure is the model's, so the camera should find it.",
  "weakly-supported",    "lower",
  "Recorded here only, neither predicted nor seen in any replicate. The thinnest positive evidence in the taxonomy, and a single record that may be a false detection, so an independent method should mostly fail to confirm it.",
  "locally absent",      "lower",
  "Not predicted, recorded in replicates, not here. Feasible in the system but not realised at this site, so the camera should not find it here.",
  "possibly forbidden",  "lower",
  "Not predicted and recorded nowhere. Either genuinely impossible or invisible to every method in use, and either way the camera should rarely find it."
)

# One category, one run: three nested models, odds ratio per SD of the posterior.
cam_fit_one <- function(d, k) {
  d$P <- d[[k]]
  na4 <- tibble(OR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_,
                OR_adj = NA_real_, p_adj = NA_real_,
                OR_adj2 = NA_real_, p_adj2 = NA_real_)
  if (dplyr::n_distinct(round(d$P, 9)) < 2) return(na4)
  d$Ps <- as.numeric(scale(d$P))
  grab <- function(form) {
    m  <- glm(as.formula(form), family = binomial, data = d)
    co <- summary(m)$coefficients
    if (!"Ps" %in% rownames(co)) return(rep(NA_real_, 4))
    c(exp(co["Ps", 1]),
      exp(co["Ps", 1] - 1.96 * co["Ps", 2]),
      exp(co["Ps", 1] + 1.96 * co["Ps", 2]),
      co["Ps", 4])
  }
  a  <- grab("camera ~ Ps")
  b  <- grab("camera ~ Ps + O_l")
  cc <- grab("camera ~ Ps + O_l + focal_site")

  # The same effect on the probability scale, which is what a reader without a
  # regression background can act on. An odds ratio is constant by
  # construction, but the probability change it implies is not: the same odds
  # ratio moves the chance less where it starts low (locally unobserved links)
  # than where it starts higher. So rather than quote one stratum, predict every
  # link twice,
  # once at the average posterior and once a standard deviation above it, and
  # average over the links. That respects the actual mix of O_l in the data and
  # answers "for a typical link here, how much does the chance move?".
  mb <- glm(camera ~ Ps + O_l, family = binomial, data = d)
  d0 <- d; d0$Ps <- 0                       # at the average posterior
  d1 <- d; d1$Ps <- 1                       # one SD above it
  pr_avg <- mean(predict(mb, d0, type = "response"))
  pr_1sd <- mean(predict(mb, d1, type = "response"))

  tibble(OR = a[1], lo = a[2], hi = a[3], p = a[4],
         OR_adj = b[1], p_adj = b[4], OR_adj2 = cc[1], p_adj2 = cc[4],
         pr_avg = pr_avg, pr_1sd = pr_1sd, pr_diff = pr_1sd - pr_avg)
}

camera_validation <- function(run_name, results_df = results, alpha = 0.05) {
  d <- results_df %>% filter(run == !!run_name, !is.na(camera))
  map_dfr(cam_expectation$category, function(k) {
    e <- cam_expectation %>% filter(category == k)
    bind_cols(
      tibble(category = k, expect = e$expect,
             links_in_category     = sum(d$det_category == k),
             confirmed_in_category = sum(d$det_category == k & d$camera == 1)),
      cam_fit_one(d, k)
    )
  }) %>%
    mutate(
      # the O_l-adjusted model decides: the posterior must beat the observation
      verdict = case_when(
        is.na(p_adj)                    ~ "not estimable",
        p_adj >= alpha                  ~ "no signal",
        expect == "unclear"             ~ "signal, no direction predicted",
        expect == "higher" & OR_adj > 1 ~ "supported",
        expect == "lower"  & OR_adj < 1 ~ "supported",
        TRUE                            ~ "CONTRADICTED"
      ),
      run = run_name, .before = 1
    )
}

cat("\n\n=====================================================================\n")
cat("CAMERA VALIDATION OF THE CATEGORIES\n")
cat("=====================================================================\n")
cam_cov <- results %>% filter(run == names(run_settings)[1], !is.na(camera))
n_run   <- sum(results$run == names(run_settings)[1])
cat(sprintf("Links inside the camera grid: %d of %d (%.0f%%); %d were recorded.\n",
            nrow(cam_cov), n_run, 100 * nrow(cam_cov) / n_run,
            sum(cam_cov$camera == 1)))

cam_report <- map_dfr(names(run_settings), camera_validation)

# Printed as one table. Column order matters here: the verdict is decided by
# the O_l-adjusted model, so OR_adj sits next to it. The raw OR is kept for
# reference but must not be read against the verdict, because it can carry the
# opposite sign, as possibly missing has done, and only the adjusted figure is
# a claim about anything the direct observation had not
# already said. n and confirmed describe the deterministic category, not the
# regression, which uses all links inside the camera grid.
cam_report %>%
  transmute(
    run, category, expect,
    n         = links_in_category,
    confirmed = confirmed_in_category,
    OR_adj    = round(OR_adj, 2),          # decisive: O_l held fixed
    p_adj     = signif(p_adj, 2),
    # the same effect as a camera-detection chance, averaged over the links
    pct_avg   = round(100 * pr_avg, 1),    # at the average posterior
    pct_1sd   = round(100 * pr_1sd, 1),    # one SD above it
    pct_diff  = round(100 * pr_diff, 1),   # percentage points
    verdict,
    OR_site   = round(OR_adj2, 2),         # also adjusted for camera placement
    p_site    = signif(p_adj2, 2),
    OR_raw    = round(OR, 2),              # reference only, see note above
    p_raw     = signif(p, 2)
  ) %>%
  arrange(match(category, cats$category), match(run, names(run_settings))) %>%
  print(n = Inf, width = Inf)

# The table above is the full record: 32 rows, one per category per run. This
# collapses it to one row per category, which is the version to read first and
# the one to lift into the SI. Two things make the collapse possible. n and
# confirmed describe the deterministic labelling, not the fit, so they are
# identical across runs and are stated once. And each run's contribution to
# the argument is a single number, OR_adj, since that is what the verdict uses.
# The question the summary answers is therefore "does this category behave as
# it claims, and do the four priors agree?", which the wide table can only be
# read down a column to answer.

SIG    <- 0.05
# The verdict counts agreement among the four PRIOR runs. The likelihood
# sensitivity run keeps its column, so its verdicts can be read off directly,
# but it does not count toward the tally: it is a near-copy of A-uniform and
# would inflate the appearance of independent replication.
N_RUNS <- length(PRIOR_RUNS)

stars <- function(p) case_when(is.na(p) ~ "",
                               p < 0.001 ~ "***",
                               p < 0.01  ~ "**",
                               p < SIG   ~ "*",
                               TRUE      ~ "")

# a range across runs, collapsed to one value when the rounded ends agree
rng_txt <- function(x, fmt = "%.0f") {
  a <- sprintf(fmt, min(x)); b <- sprintf(fmt, max(x))
  if (a == b) a else paste0(a, "-", b)
}

cam_by_run <- cam_report %>%
  filter(run %in% PRIOR_RUNS) %>%
  mutate(cell = ifelse(is.na(OR_adj), "--",
                       sprintf("%.2f%s", OR_adj, stars(p_adj)))) %>%
  select(category, run, cell) %>%
  pivot_wider(names_from = run, values_from = cell)

# Eight odds-ratio columns would be unreadable, and the question the twins
# answer is narrow: does removing nu change the CALL? So each twin is compared
# with its own base run and the four comparisons collapse to one column. A call
# is the pair (significant or not, and which way), which is exactly what the
# verdict is built from.
cam_no_nu <- cam_report %>%
  mutate(call = case_when(is.na(p_adj) ~ "na",
                          p_adj >= SIG ~ "none",
                          OR_adj > 1   ~ "up",
                          TRUE         ~ "down"),
         base = sub(" \\(no nu\\)$", "", run)) %>%
  # exactly two rows per (category, base): a run and its twin. They agree when
  # the two calls are the same, so no reshaping is needed.
  group_by(category, base) %>%
  summarise(same = n_distinct(call) == 1, .groups = "drop_last") %>%
  summarise(no_nu = sprintf("%d/%d same", sum(same), n()), .groups = "drop")

cam_summary <- cam_report %>%
  mutate(
    # significant AND running the way the category claims, or against it
    right_way = p_adj < SIG & ((expect == "higher" & OR_adj > 1) |
                               (expect == "lower"  & OR_adj < 1)),
    wrong_way = p_adj < SIG & ((expect == "higher" & OR_adj < 1) |
                               (expect == "lower"  & OR_adj > 1))
  ) %>%
  group_by(category, expect) %>%
  summarise(
    n         = first(links_in_category),
    confirmed = first(confirmed_in_category),
    pct_conf  = round(100 * first(confirmed_in_category) /
                            first(links_in_category)),
    n_sig     = sum(p_adj < SIG & run %in% PRIOR_RUNS,   na.rm = TRUE),
    supports  = sum(right_way    & run %in% PRIOR_RUNS,  na.rm = TRUE),
    against   = sum(wrong_way    & run %in% PRIOR_RUNS,  na.rm = TRUE),
    # of the runs that found a directional effect, how many keep it once
    # camera placement is adjusted for as well
    site_ok   = sum((right_way | wrong_way) & p_adj2 < SIG &
                    run %in% PRIOR_RUNS, na.rm = TRUE),
    # the probability translation, ranged over the runs that found an effect.
    # scalar if(), not ifelse(), because ifelse evaluates both branches and
    # the range of an empty selection would be -Inf
    chance    = local({
      keep <- !is.na(p_adj) & p_adj < SIG & run %in% PRIOR_RUNS
      if (!any(keep)) "--" else
        sprintf("%s -> %s", rng_txt(100 * pr_avg[keep]),
                            rng_txt(100 * pr_1sd[keep]))
    }),
    .groups   = "drop"
  ) %>%
  mutate(
    verdict = case_when(
      expect == "unclear"         ~ sprintf("no direction predicted (%d/%d signal)",
                                            n_sig, N_RUNS),
      against > 0 & supports == 0 ~ sprintf("CONTRADICTED in %d/%d", against, N_RUNS),
      against > 0                 ~ sprintf("mixed: %d for, %d against",
                                            supports, against),
      supports > 0                ~ sprintf("supported in %d/%d", supports, N_RUNS),
      TRUE                        ~ sprintf("no signal (0/%d)", N_RUNS)),
    site = case_when(
      supports + against == 0        ~ "",
      site_ok == supports + against  ~ "holds",
      TRUE ~ sprintf("%d of %d", site_ok, supports + against))
  )

cam_summary_tbl <- cam_summary %>%
  left_join(cam_by_run, by = "category") %>%
  left_join(cam_no_nu,  by = "category") %>%
  select(category, expect, n, confirmed, pct_conf,
         all_of(PRIOR_RUNS), verdict, site, no_nu, chance) %>%
  arrange(match(category, cats$category))

# COLUMNS
#   n, confirmed  size of the DETERMINISTIC category among links inside the
#                 camera grid, and how many of those the cameras recorded, with
#                 pct_conf the share. Not the size of the whole category.
#                 Identical in every run, so stated once. Descriptive only: the
#                 regressions use all links inside the camera grid.
#   run columns   OR_adj, the odds ratio per 1 SD of that category's posterior
#                 with the direct observation held fixed. Above 1, a higher
#                 posterior goes with MORE camera records.
#                 *** p<0.001   ** p<0.01   * p<0.05
#   verdict       how many runs move the way the category claims. Read it with
#                 expect, which is the ecological claim, not the z_l bit.
#   site          of those runs, how many survive adjusting for camera
#                 placement too. Blank where no run found a direction.
#   no_nu         how many of the four priors give the SAME call once link
#                 heterogeneity is removed. A call is "significant and upward",
#                 "significant and downward", or "no signal", which is what the
#                 verdict is built from. 4/4 means nu changes nothing here.
#   chance        camera-record chance in percent at the average posterior ->
#                 one SD above it, ranged over the significant runs. This is
#                 the figure to quote in prose, not the odds ratio.
cat("\n\nSUMMARY: one row per category, the four runs side by side\n\n")
print(cam_summary_tbl, n = Inf, width = Inf)

write_csv(cam_report, file.path(OUT_DIR, "camera_validation_by_category.csv"))
write_csv(cam_summary_tbl, file.path(OUT_DIR, "camera_validation_summary.csv"))
cat(sprintf("\nWritten: %s\n         %s\n",
            file.path(OUT_DIR, "camera_validation_by_category.csv"),
            file.path(OUT_DIR, "camera_validation_summary.csv")))


# ---- 10. Evidence accumulation on the empirical rates (Fig. 3a applied) ----
#
# Figure 3a of the manuscript shows confidence climbing as replicates keep
# recording a link, and stopping at kappa. That figure uses illustrative rates.
# This is the same construction driven by the rates fitted here, so it says
# what replication is actually worth in this system.
#
# CONSTRUCTION
# The evidence is the possibly-missing signature: predicted (Y = 1), not seen
# locally (O_l = 0), and detected in every replicate that could have found it
# (n = R). R is then swept. Curves come from posterior() itself rather than a
# second implementation, so this figure and the analysis cannot drift apart.
#
# WHAT THE SOLID AND DASHED PARTS MEAN
#   solid   R <= 5, which six sites can deliver. This is the study as built.
#   dashed  R > 5, hypothetical extra replicate sites that do not exist.
# The dashed section is the answer to "how many more sites would we need?": the
# closer the solid line already sits to kappa at the design limit, the less
# extra sites can add.
#
# WHAT THE FIGURE SHOWS ON THESE RATES
# How fast the curve climbs is set by how far p1 sits above f: the wider the
# gap, the stronger each recording replicate is as evidence. The value it
# climbs to is kappa, and kappa does not involve p1 at all. The ceiling is a
# property of the model and the local method, not of the sampling design, so
# the way past it is a better model (lower eps_Y) or a better local method
# (lower eps_l), which is exactly the action the taxonomy prescribes. The
# printout below gives the posterior at the design limit, to set against kappa.
#
# COLOURS
# Each category keeps the colour it carries in the maps, which frees linetype
# to mean observed versus hypothetical. The manuscript figure uses linetype to
# pair categories instead, because it has no extrapolation to mark.

R_MAX_DESIGN <- N_REPLICATES   # every other site, the most the design provides
R_SHOW       <- 12     # far enough past saturation to make the plateau obvious

accum <- map_dfr(0:R_SHOW, function(R) {
  p <- posterior(Y = 1, O_l = 0, n = R, R = R)      # detected in every replicate
  tibble(R = R, category = colnames(p), P = as.numeric(p))
})

KAPPA_EMP <- kappa(EPS_Y, EPS_L, F_POS)

# the four categories of Fig. 3a: the two the evidence is deciding between,
# and the two they would become if the local axis had erred
FIG3A_CATS <- c("possibly missing", "phantom", "recurrent", "locally unique")

acc_plot <- accum %>%
  filter(category %in% FIG3A_CATS) %>%
  mutate(category = factor(category, levels = FIG3A_CATS))

# feasibility confidence: realised here, or realisable in the replicates
phi_emp <- accum %>%
  group_by(R) %>%
  summarise(phi = 1 - sum(P[category %in% NOT_FEASIBLE]), .groups = "drop")

# the replicate count at which each category first reaches a target
r_needed <- function(cat, target) {
  s <- acc_plot %>% filter(category == cat, P >= target)
  if (nrow(s) == 0) NA_integer_ else min(s$R)
}
TARGET <- 0.50
R_STAR <- r_needed("possibly missing", TARGET)

cat("\n\nEVIDENCE ACCUMULATION ON THE FITTED RATES\n")
cat(sprintf("  kappa = %.3f, the ceiling for a link not recorded locally\n",
            KAPPA_EMP))
cat(sprintf("  P(possibly missing) reaches %.2f at R = %s, and %.3f at the design limit R = %d\n",
            TARGET, ifelse(is.na(R_STAR), "never", R_STAR),
            acc_plot$P[acc_plot$category == "possibly missing" &
                         acc_plot$R == R_MAX_DESIGN],
            R_MAX_DESIGN))
print(acc_plot %>% filter(R <= 6) %>%
        pivot_wider(names_from = category, values_from = P) %>%
        mutate(across(-R, ~round(.x, 3))) %>% as.data.frame(), row.names = FALSE)

si_base10 <- theme_classic(base_size = 10) +
  theme(axis.line       = element_line(colour = "grey40", linewidth = 0.3),
        axis.ticks      = element_line(colour = "grey40", linewidth = 0.3),
        legend.position = "bottom",
        legend.title    = element_blank(),
        legend.key.width = unit(18, "pt"),
        legend.margin   = margin(t = 0, b = 0))

fig_accum <- ggplot(acc_plot, aes(R, P, colour = category)) +
  # the ceiling, and the design limit
  geom_hline(yintercept = KAPPA_EMP, linetype = "22", colour = "grey45",
             linewidth = 0.4) +
  annotate("text", x = R_SHOW, y = KAPPA_EMP, label = "kappa == 0.533",
           parse = TRUE, hjust = 1, vjust = -0.6, size = 3, colour = "grey35",
           fontface = "bold") +
  geom_vline(xintercept = R_MAX_DESIGN, linetype = "dotted",
             colour = "grey60", linewidth = 0.35) +
  annotate("text", x = R_MAX_DESIGN, y = 0.02, label = "six sites",
           hjust = -0.1, size = 2.8, colour = "grey45") +
  # feasibility confidence, as in the manuscript figure
  geom_line(data = phi_emp, aes(R, phi), inherit.aes = FALSE,
            colour = "#2A9D8F", linewidth = 0.8, linetype = "12") +
  annotate("text", x = R_SHOW, y = 1, label = "feasibility confidence",
           hjust = 1, vjust = 1.6, size = 2.9, colour = "#2A9D8F",
           fontface = "bold") +
  # solid where the design reaches, dashed where it does not
  geom_line(data = ~ filter(.x, R <= R_MAX_DESIGN), linewidth = 0.85) +
  geom_line(data = ~ filter(.x, R >= R_MAX_DESIGN), linewidth = 0.85,
            linetype = "22") +
  scale_colour_manual(values = CATEGORY_COLOUR[FIG3A_CATS],
                      labels = c("possibly missing (1,0,1)", "phantom (1,0,0)",
                                 "recurrent (1,1,1)", "locally unique (1,1,0)")) +
  scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(round(100*x), "%")) +
  scale_x_continuous(breaks = seq(0, R_SHOW, 2)) +
  labs(x = "Replicates recording the link, R",
       y = "Posterior probability",
       # two lines, because one would run past the panel and be clipped
       caption = paste0(
         "Solid: within the six-site design. Dashed: hypothetical extra replicates.\n",
         # plain ASCII: the pdf device cannot encode Greek or subscripts here
         # and silently substitutes dots for them
         "Fitted rates eps_Y = ", round(EPS_Y, 2),
         ", eps_l = ", round(EPS_L, 2),
         ", f = ", F_POS, ", p1 = ", round(P1, 2),
         ", nu = ", round(NU, 2), ".")) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  si_base10 +
  theme(plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey35",
                                    lineheight = 1.2))

ggsave(file.path(OUT_DIR, "fig_accumulation_empirical.pdf"), fig_accum,
       width = 6.5, height = 5)
ggsave(file.path(OUT_DIR, "fig_accumulation_empirical.png"), fig_accum,
       width = 6.5, height = 5, dpi = 300)
if (interactive()) print(fig_accum)
cat(sprintf("\nWritten: %s\n",
            file.path(OUT_DIR, "fig_accumulation_empirical.pdf")))


# ---- 11. A per-species view ----
# PURPOSE   Field effort is often planned per species. How many of each plant's
#           candidate links does the posterior place in each category?
# WHAT      Per plant at the richest site, sum each category's posterior over
#           the plant's links (expected count) and give its SD.
# DECISIONS The sum, not the mean share or a count of most-likely categories.
#           B-uniform and B-degree, since the prior moves the posterior far more
#           than where the model enters. nu kept. Rationale and caveats are in
#           the note, section "Where to look, species by species".

SPECIES_RUNS <- c("B-uniform", "B-degree")

per_plant <- results %>%
  filter(run %in% SPECIES_RUNS, focal_site == richest) %>%
  select(run, focal_site, plant, det_category, all_of(CAT_ORDER)) %>%
  pivot_longer(all_of(CAT_ORDER), names_to = "category", values_to = "p") %>%
  group_by(run, focal_site, plant, category) %>%
  summarise(candidates    = n(),
            deterministic = sum(det_category == category),
            expected      = sum(p),
            sd            = sqrt(sum(p * (1 - p))),   # links taken as independent
            share         = expected / candidates,     # average posterior
            .groups = "drop")

write_csv(per_plant,file.path(OUT_DIR, "per_plant_expected_counts.csv"))

# a plant's eight expected counts must add up to its candidate links
sum_gap <- per_plant %>%
  group_by(run, plant) %>%
  summarise(gap = abs(sum(expected) - first(candidates)), .groups = "drop")
cat(sprintf("  [%s] per plant, expected counts sum to the candidate links (max gap %.1e)\n",
            ifelse(max(sum_gap$gap) < 1e-8, "OK", "!!"), max(sum_gap$gap)))

# same ordering rule as make_map: this site's plants by pooled degree
pla_site <- per_plant %>% distinct(plant) %>%
  left_join(deg_pla, by = c("plant" = "lower_level")) %>%
  arrange(desc(deg)) %>% pull(plant)

species_figure <- function(run_name) {
  d <- per_plant %>%
    filter(run == run_name) %>%
    mutate(plant    = factor(plant, levels = pla_site),
           category = factor(category, levels = CAT_ORDER))

  ggplot(d, aes(y = plant)) +
    # whisker is expected +/- SD, cut at zero
    geom_segment(aes(x = pmax(expected - sd, 0), xend = expected + sd,
                     yend = plant), colour = "grey45", linewidth = 0.5) +
    geom_point(aes(x = deterministic), shape = 124, size = 3,
               colour = "grey15") +
    geom_point(aes(x = expected, fill = category), shape = 21, size = 2.4,
               colour = "grey25", stroke = 0.3) +
    # predicted four on the top row, unpredicted four below
    facet_wrap(~ category, nrow = 2, scales = "free_x") +
    scale_fill_manual(values = CATEGORY_COLOUR) +
    scale_y_discrete(labels = ital) +
    scale_x_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0.02, 0.08))) +
    # plain ASCII: the pdf device cannot encode the plus-minus sign
    labs(title = sprintf("%s, %s", richest, run_name),
         x = "Links per plant: expected +/- SD (dot and whisker), deterministic count (|)",
         y = NULL) +
    si_base +
    theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
          plot.title         = element_text(size = 10, face = "bold"))
}

for (rn in SPECIES_RUNS) {
  fig_sp <- species_figure(rn)
  stem <- file.path(OUT_DIR, paste0("fig_per_plant_", run_stem(rn)))
  ggsave(paste0(stem, ".pdf"), fig_sp, width = 12, height = 7)
  ggsave(paste0(stem, ".png"), fig_sp, width = 12, height = 7, dpi = 300)
  if (interactive()) print(fig_sp)
}

cat("\nEXPECTED LINKS PER PLANT (richest site)\n")
per_plant %>%
  mutate(expected = round(expected, 1),
         category = factor(category, levels = CAT_ORDER)) %>%
  arrange(category) %>%
  select(run, plant, category, expected) %>%
  pivot_wider(names_from = category, values_from = expected) %>%
  arrange(run, match(plant, pla_site)) %>%
  print(n = Inf, width = Inf)

# --- step 1, heatmap: each plant's average posterior ------------------------
# PURPOSE   See every plant's composition at once.
# WHAT      Plants x categories, cell = average posterior over the plant's
#           links (expected count / candidates), so each row sums to 1.
# DECISIONS Each column runs white -> its category colour from the map key, on
#           one intensity scale for all cells, so a row is read across
#           categories. Values printed, as small shares are hard to read by hue
#           and pale category colours never get dark.

share_gap <- per_plant %>% group_by(run, plant) %>%
  summarise(gap = abs(sum(share) - 1), .groups = "drop")
cat(sprintf("  [%s] per plant, average posteriors sum to 1 (max gap %.1e)\n",
            ifelse(max(share_gap$gap) < 1e-8, "OK", "!!"), max(share_gap$gap)))

share_heatmap <- function(run_name) {
  top <- max(per_plant$share[per_plant$run == run_name])

  d <- per_plant %>%
    filter(run == run_name) %>%
    mutate(plant    = factor(plant, levels = pla_site),
           category = factor(category, levels = CAT_ORDER),
           # white -> category colour, as in make_map, one intensity scale
           fill = map2_chr(as.character(category), share / top, function(k, v)
             colorRampPalette(c("white", CATEGORY_COLOUR[[k]]))(101)[round(100 * v) + 1]),
           # dark text on light cells, white on dark, by perceived luminance
           lum  = colSums(col2rgb(fill) * c(0.299, 0.587, 0.114)) / 255,
           ink  = if_else(lum < 0.5, "white", "grey15"))

  ggplot(d, aes(category, plant)) +
    geom_tile(aes(fill = fill), colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", share), colour = ink), size = 2.8) +
    scale_fill_identity() +
    scale_colour_identity() +
    scale_x_discrete(position = "top") +
    scale_y_discrete(labels = ital) +
    labs(title = sprintf("%s, %s", richest, run_name), x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 10) +
    theme(panel.grid  = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 0, vjust = 0),
          plot.title  = element_text(size = 10, face = "bold"),
          # room for the last angled column label
          plot.margin = margin(5, 60, 5, 5))
}

for (rn in SPECIES_RUNS) {
  fig_hm <- share_heatmap(rn)
  stem <- file.path(OUT_DIR, paste0("fig_per_plant_share_", run_stem(rn)))
  ggsave(paste0(stem, ".pdf"), fig_hm, width = 8, height = 5)
  ggsave(paste0(stem, ".png"), fig_hm, width = 8, height = 5, dpi = 300)
  if (interactive()) print(fig_hm)
}

# --- step 2: one plant, where its posterior mass goes ------------------------
# PURPOSE   For a chosen plant, compare its categories with each other.
# WHAT      The step-1 values for that plant, one row per category, one figure
#           per run.
# DECISIONS One axis for all eight rows, and the same axis in both runs'
#           figures, so categories and runs compare directly. Rationale in the
#           note.

# choose the plants to draw; any plant at the richest site works
FOCAL_PLANTS <- c("Sedum sediforme", "Teucrium capitatum")

stopifnot(all(FOCAL_PLANTS %in% pla_site))

plant_figure <- function(plant_name, run_name) {
  d_all <- per_plant %>% filter(plant == plant_name)
  # shared by both runs, so their figures sit on one scale
  x_max <- max(d_all$expected + d_all$sd, d_all$deterministic)

  d <- d_all %>%
    filter(run == run_name) %>%
    mutate(category = factor(category, levels = rev(CAT_ORDER)))

  ggplot(d, aes(y = category)) +
    geom_segment(aes(x = pmax(expected - sd, 0), xend = expected + sd,
                     yend = category), colour = "grey45", linewidth = 0.5) +
    geom_point(aes(x = deterministic), shape = 124, size = 3.2,
               colour = "grey15") +
    geom_point(aes(x = expected, fill = category), shape = 21, size = 2.8,
               colour = "grey25", stroke = 0.3) +
    scale_fill_manual(values = CATEGORY_COLOUR) +
    scale_x_continuous(limits = c(0, x_max),
                       expand = expansion(mult = c(0.02, 0.06))) +
    labs(title = bquote(italic(.(plant_name)) ~ "at" ~ .(richest) * "," ~
                          .(run_name) * "," ~ .(unique(d$candidates)) ~
                          "candidate links"),
         x = "Links: expected +/- SD (dot and whisker), deterministic count (|)",
         y = NULL) +
    si_base +
    theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
          plot.title         = element_text(size = 10, face = "bold"))
}

for (pl in FOCAL_PLANTS) for (rn in SPECIES_RUNS) {
  fig_pl <- plant_figure(pl, rn)
  stem <- file.path(OUT_DIR, sprintf("fig_plant_%s_%s", run_stem(pl), run_stem(rn)))
  ggsave(paste0(stem, ".pdf"), fig_pl, width = 6.5, height = 4)
  ggsave(paste0(stem, ".png"), fig_pl, width = 6.5, height = 4, dpi = 300)
  if (interactive()) print(fig_pl)
}


# ---- 12. Category leakage ----
# PURPOSE   Where does the posterior move each deterministic category once
#           uncertainty is incorporated?
# WHAT      A network of the eight categories. The arrow k -> m carries
#           F(k -> m) = sum over links assigned k of P(m | E), the expected
#           number of links moving; what stays is F(k -> k). All sites.
# DECISIONS Expected links, not counts of changed most-likely categories. All
#           arrows, absolute or relative weights (share of the source category),
#           one width scale per weighting for every run. Layout of
#           the manuscript's taxonomy figure: columns = observed elsewhere,
#           rows = TP, FP | TN, FN, so a move to a neighbour doubts one bit.
#           Arrow colour = source category. Rationale in the note.

library(ggraph)
library(tidygraph)

LEAK_RUNS <- SPECIES_RUNS

# node positions: the taxonomy figure, predicted block above unpredicted
LEAK_W <- 11
leak_nodes <- tibble(
  category = c("locally unique", "recurrent", "phantom", "possibly missing",
               "possibly forbidden", "locally absent", "weakly-supported",
               "model-elusive"),
  row = c(1, 1, 2, 2, 3, 3, 4, 4)) %>%
  left_join(cats, by = "category") %>%
  mutate(x = zr * LEAK_W,
         y = c(9, 6.2, 2.8, 0)[row],
         signature = sprintf("(%d,%d,%d)", zY, zl, zr))

# F(k -> m) for one run, including what stays (k == m)
leak_flows <- function(run_name) {
  results %>%
    filter(run == run_name) %>%
    select(det_category, all_of(CAT_ORDER)) %>%
    pivot_longer(all_of(CAT_ORDER), names_to = "to", values_to = "p") %>%
    group_by(from = det_category, to) %>%
    summarise(w = sum(p), .groups = "drop")
}

# share = F(k -> m) / N_k^det, the share of category k's links moving to m
flows_all <- map_dfr(LEAK_RUNS, ~ leak_flows(.x) %>% mutate(run = .x)) %>%
  group_by(run, from) %>%
  mutate(share = w / sum(w)) %>%
  ungroup()
write_csv(flows_all, file.path(OUT_DIR, "category_leakage_flows.csv"))

# weights  "absolute": expected links moved, which the largest categories
#                       dominate simply by size
#          "relative": share of the source category's links moved, which shows
#                       how stable each category is regardless of size
# The largest moving flow in any run sets each shared width scale.
LEAK_SCALE <- list(
  absolute = list(col = "w",     max = max(flows_all$w[flows_all$from != flows_all$to]),
                  key = c(5, 25, 100), key_lab = c("5", "25", "100"),
                  title = "expected links moved"),
  relative = list(col = "share", max = max(flows_all$share[flows_all$from != flows_all$to]),
                  key = c(0.05, 0.15, 0.3), key_lab = c("5%", "15%", "30%"),
                  title = "share of the category's links moved"))

leakage_network <- function(run_name, weights = c("absolute", "relative")) {
  weights <- match.arg(weights)
  sc <- LEAK_SCALE[[weights]]
  flows <- filter(flows_all, run == run_name) %>%
    mutate(weight = .data[[sc$col]])

  nodes <- leak_nodes %>%
    left_join(flows %>% group_by(category = from) %>% summarise(det = sum(w)),
              by = "category") %>%
    left_join(flows %>% group_by(category = to) %>% summarise(expected = sum(w)),
              by = "category") %>%
    left_join(flows %>% filter(from == to) %>% select(category = from, kept = w),
              by = "category") %>%
    mutate(across(c(det, expected, kept), ~ replace_na(.x, 0)),
           diam   = 5 + 15 * sqrt(det / max(det)),        # mm, area ~ count
           detail = sprintf("%s  ·  %s → %s  ·  %s%% kept",
                            signature, format(round(det), big.mark = ","),
                            format(round(expected), big.mark = ","),
                            round(100 * kept / pmax(det, 1))))

  edges <- flows %>%
    filter(from != to) %>%
    left_join(select(leak_nodes, from = category, x0 = x, y0 = y, r0 = row), by = "from") %>%
    left_join(select(leak_nodes, to = category, x1 = x, y1 = y, r1 = row), by = "to") %>%
    mutate(
      # a move within a column that skips a row bends inward, clear of the
      # nodes it passes; the rest bend slightly left, which separates the two
      # directions of a pair. Positive strength bends left of travel.
      long     = x0 == x1 & abs(r1 - r0) >= 2,
      out_x    = ifelse(x0 == 0, -1, 1),
      down     = sign(y1 - y0),
      strength = round(ifelse(long, out_x * down * ifelse(down < 0, 0.30, 0.42), 0.14), 2),
      src      = from,
      from     = match(src, nodes$category),
      to       = match(to, nodes$category),
      cap_s    = nodes$diam[from],
      cap_e    = nodes$diam[to]) %>%
    arrange(weight)                                    # heavy arrows on top

  g <- tbl_graph(nodes = nodes,
                 edges = select(edges, from, to, weight, src, strength, cap_s, cap_e))

  p <- ggraph(g, layout = "manual", x = nodes$x, y = nodes$y)
  # strength is a layer parameter, so one layer per value; !! fixes s now,
  # otherwise every layer would filter on the loop's last value
  for (s in sort(unique(edges$strength))) {
    p <- p + geom_edge_arc(
      aes(width = weight, colour = src, filter = strength == !!s,
          start_cap = circle(cap_s / 2 + 0.8, "mm"),
          end_cap   = circle(cap_e / 2 + 1.2, "mm")),
      strength = s, alpha = 0.6, lineend = "butt",
      arrow = arrow(length = unit(3, "mm"), type = "closed", angle = 25))
  }

  width_range <- c(0.15, 5)
  key <- tibble(weight = sc$key, label = sc$key_lab,
                x = LEAK_W / 2 - 2.2 + (0:2) * 2.2, y = -2.6)

  # labels on the outer side of each column, clear of the arrows
  MM_PER_UNIT <- 10.5                   # approximate, at the saved size below
  lab <- nodes %>%
    mutate(side = ifelse(zr == 0, -1, 1),
           lx   = x + side * (diam / 2 / MM_PER_UNIT + 0.45),
           hj   = ifelse(side < 0, 1, 0))

  p +
    geom_node_point(aes(size = I(diam), fill = category), shape = 21,
                    colour = "white", stroke = 1.4) +
    geom_text(data = lab, aes(x = lx, y = y + 0.12, label = category, hjust = hj),
              vjust = 0, size = 3.6, fontface = "bold", colour = "#2b2b3a",
              inherit.aes = FALSE) +
    geom_text(data = lab, aes(x = lx, y = y - 0.12, label = detail, hjust = hj),
              vjust = 1, size = 2.6, colour = "#707080", inherit.aes = FALSE) +
    annotate("text", x = c(0, LEAK_W), y = 10.6,
             label = c("not observed elsewhere", "observed elsewhere"),
             size = 3.2, colour = "#8a8a9a") +
    annotate("text", x = -7.1, y = c(7.6, 1.4),
             label = c("predicted", "not predicted"),
             size = 3.2, colour = "#8a8a9a", fontface = "italic", angle = 90) +
    annotate("segment", x = -6.75, xend = -6.75, y = c(5.9, -0.3),
             yend = c(9.3, 3.1), colour = "#d0d0de", linewidth = 0.5) +
    # width key, drawn on the same scale as the arrows
    geom_segment(data = key,
                 aes(x = x - 0.6, xend = x + 0.4, y = y, yend = y,
                     linewidth = I(scales::rescale(weight, width_range, c(0, sc$max)) * 0.75)),
                 colour = "#9a9aae", lineend = "butt", inherit.aes = FALSE) +
    geom_text(data = key, aes(x = x + 0.6, y = y, label = label), hjust = 0,
              size = 2.8, colour = "#707080", inherit.aes = FALSE) +
    annotate("text", x = LEAK_W / 2, y = -1.9, label = sc$title,
             size = 2.8, colour = "#707080") +
    scale_edge_colour_manual(values = CATEGORY_COLOUR, guide = "none") +
    scale_fill_manual(values = CATEGORY_COLOUR, guide = "none") +
    scale_edge_width(range = width_range, limits = c(0, sc$max), guide = "none") +
    coord_fixed(clip = "off", xlim = c(-7.4, LEAK_W + 5.2), ylim = c(-2.9, 10.9)) +
    labs(title = run_name) +
    theme_void(base_size = 11) +
    theme(plot.title  = element_text(face = "bold", colour = "#2b2b3a", hjust = 0.02),
          plot.margin = margin(10, 20, 10, 20))
}

for (rn in LEAK_RUNS) for (wt in c("absolute", "relative")) {
  fig_leak <- leakage_network(rn, weights = wt)
  stem <- file.path(OUT_DIR, sprintf("fig_category_leakage_%s%s", run_stem(rn),
                                     if (wt == "relative") "_relative" else ""))
  # cairo keeps the arrow and dot glyphs the default pdf device drops
  ggsave(paste0(stem, ".pdf"), fig_leak, width = 10.5, height = 6.9,
         device = cairo_pdf)
  ggsave(paste0(stem, ".png"), fig_leak, width = 10.5, height = 6.9, dpi = 300,
         bg = "white")
  if (interactive()) print(fig_leak)
}

cat("\nCATEGORY LEAKAGE: share of each deterministic category kept\n")
flows_all %>%
  group_by(run, category = from) %>%
  summarise(kept = sum(w[from == to]) / sum(w), .groups = "drop") %>%
  mutate(kept = round(kept, 2)) %>%
  pivot_wider(names_from = run, values_from = kept) %>%
  arrange(match(category, CAT_ORDER)) %>%
  print()
