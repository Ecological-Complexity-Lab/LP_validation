# =============================================================================
# Bayesian categorisation of predicted links: the Cabrera worked example
#
# Every candidate link gets eight posterior probabilities, one per category,
# instead of a single hard label.
#
#   1. Load the data and assemble the evidence
#   2. Estimate the rates
#   3. Calibrate the model score, and check it against the hard predictions
#   4. The posterior function
#   5. Check the function against every number published in the SI
#   6. Run it
#   7. The map
#
# There is one posterior function. Every SI section is a setting of its
# arguments, not a separate implementation, so the equation printed in the
# manuscript and the code producing the results are the same object.
# =============================================================================

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


# =============================================================================
# 1. THE EVIDENCE
# =============================================================================
# Each link needs four things: Y, O_l, n and R.
#
# On R. Every site is a replicate, so R = 5 everywhere. Co-occurrence lives in
# the realisation rate instead, as SI Table S6 prescribes: rho(k) = c(k) * r,
# where c(k) = 1 only where both partners are present. At a site where they do
# not co-occur, p1 and p0 are both zero, that replicate contributes the same
# factor of 1 to all eight categories, and it cancels. So the arithmetic
# reduces to counting only the sites where the pair could have been found,
# which is what R below holds. 582 of 1624 links have R = 0: their replicate
# evidence is silent rather than negative, and the posterior will say so.

raw <- read_csv(DATA_OBS, show_col_types = FALSE)

obs <- raw %>%
  filter(method == "obs") %>%
  mutate(pair = paste(higher_level, lower_level, sep = "||"))

# camera records, used twice: to estimate the miss rate, and to mark the map
cam <- raw %>%
  filter(method == "rpi") %>%
  transmute(higher_level, lower_level, focal_site, camera = ground_truth)

ev <- obs %>%
  group_by(pair) %>%
  mutate(
    R = n_distinct(focal_site) - 1,        # co-occurring replicate sites
    n = sum(ground_truth) - ground_truth   # detections among them
  ) %>%
  ungroup() %>%
  transmute(
    focal_site,
    pollinator = higher_level,
    plant      = lower_level,
    Y   = prediction,                      # hard model call
    O_l = ground_truth,                    # seen here?
    n, R,
    score = qlogis(probability)            # raw softImpute value, see step 3
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


# =============================================================================
# 2. THE RATES
# =============================================================================

# --- eps_Y: model error -----------------------------------------------------
# 1 - balanced accuracy. This equals the mean of the two directional error
# rates, so it is a rate measured at the threshold actually used, and it is
# symmetric, which is what SI Table S4 asks for. Table S1 defines this error
# against the observations; whether the observations are the truth is handled
# downstream by eps_l and f, so it is not also a concern here.

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
# A lower bound, because the cameras miss links too.

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
# p1 = rho * (1 - eps_l), Eq. S4. Only the product ever enters the replicate
# factor, so we estimate it directly and never assume rho (we cannot separate rho from eps_l).
#
# Fit to pairs seen at least once, which is our working definition of
# "realisable in the replicates". Conditioning on n >= 1 truncates the
# distribution, so the naive proportion is biased upward. The zero-truncated
# likelihood corrects for it. The beta-binomial version also returns nu, the
# concentration that SI Section S7 needs.

estimate_p1 <- function(obs) {
  pairs <- obs %>%
    group_by(pair) %>%
    summarise(R = n_distinct(focal_site), n = sum(ground_truth),
              .groups = "drop") %>%
    filter(R >= 2, n >= 1)

  # log P(n | R) for a beta-binomial, without the binomial coefficient (it
  # cancels everywhere), minus log P(n >= 1) to correct the truncation
  nll <- function(par) {
    m <- plogis(par[1]); nu <- exp(par[2])
    a <- m * nu; b <- (1 - m) * nu
    lp <- function(k, N) lbeta(k + a, N - k + b) - lbeta(a, b)
    -sum(lp(pairs$n, pairs$R) - log1p(-exp(lp(0, pairs$R))))
  }

  fit <- optim(c(0, log(5)), nll)
  list(p1 = plogis(fit$par[1]), nu = exp(fit$par[2]),
       n_pairs = nrow(pairs),
       n_all   = n_distinct(obs$pair),
       set     = "species pairs co-occurring at 2+ sites and recorded at 1+ of them")
}

rep_fit <- estimate_p1(obs)
P1  <- rep_fit$p1
NU  <- rep_fit$nu
P0  <- F_POS
RHO <- P1 / (1 - EPS_L)     # derived and reported, never an input

cat("\nRATES\n")
cat(sprintf("  eps_Y = %.3f   (1 - balanced accuracy)\n", EPS_Y))
cat(sprintf("  eps_l = %.3f   (camera contrast)\n", EPS_L))
cat(sprintf("  f     = %.3f   (fixed a priori)\n", F_POS))
cat(sprintf("  p1    = %.3f   nu = %.2f\n", P1, NU))
cat(sprintf("          estimated on %d of %d species pairs: %s\n",
            rep_fit$n_pairs, rep_fit$n_all, rep_fit$set))
cat("          Zero-truncated because conditioning on 'recorded at least once'\n")
cat("          selects for high p1; the truncation corrects that upward bias.\n")

# =============================================================================
# 3. CALIBRATING THE MODEL SCORE
# =============================================================================
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
# asserting it.
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
cat("       So the two scales give the SAME classification. Only the number\n")
cat(sprintf("       on the threshold changes, from 0.70 to %.3f.\n", Q_THRESHOLD))

# How far apart are the six fold-specific calibrations? Small spread means the
# score behaves the same way in rocky, dune and pine sites, so no site needs
# special handling. This is information, not a pass/fail test.
q_global <- predict(glm(O_l ~ score, family = binomial, data = ev),
                    type = "response")

cat(sprintf("  info calibration varies across sites by at most %.3f in q\n",
            max(abs(ev$q - q_global))))
cat(sprintf("       (pooled Spearman is %.4f rather than 1 because this is six\n",
            cor(ev$score, ev$q, method = "spearman")))
cat("        models, not one; that is expected and not a failure)\n")

cat(sprintf("  note %d links are predicted (Y = 1) yet calibrate below 0.5.\n",
            sum(ev$Y == 1 & ev$q < 0.5)))
cat("       These are where the threshold discards information: a precision-\n")
cat("       weighted rule flags them, a calibrated reading puts them just on\n")
cat("       the absent side. Run B treats them differently from Run A, which\n")
cat("       is the point of SI Fig. S9.\n")

# =============================================================================
# 4. THE POSTERIOR
# =============================================================================
# Calculated from the master equation the unifies all measurable elements described in the SI. 
# Nothing is hardcoded: every SI section is a setting of the arguments, not a separate branch.
#
#   Section S2  eps_r supplied, count = FALSE   symmetric binary evidence
#   Section S5  eps_r = NULL,   count = FALSE   directional rates from p1 and f
#   Section S6  count = TRUE,   nu = Inf        cumulative binomial
#   Section S7  count = TRUE,   nu finite       heterogeneous p1
#   Section S8  s = 0, pi_Y = q                 model as prior
#
# The binary forms are not an alternative formula. They are the same replicate
# distribution with the count aggregated: P(n = 0) against P(n >= 1) instead of
# P(n = k). Setting count = FALSE coarsens the evidence; it does not change the
# model. Likewise nu = Inf is the limit in which all links share one p1.
#
# eps_r is the one genuine override. Section S2 treats the replicate error as a
# free-standing rate, before it is grounded in rho and eps_l, so no value of p1
# reproduces it. Leave it NULL for everything from Section S5 onwards.
#
# Arguments
#   Y, O_l, n, R      evidence; vectors of equal length
#   eps_Y             model error; scalar, or c(false_neg, false_pos)
#   eps_l, f          local miss rate and false-detection rate
#   p1bar, nu         mean and concentration of p1 across links; nu = Inf for
#                     one shared rate
#   p0                replicate false detection; defaults to f as Table S4
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
    #     non-realisable categories keep a fixed rate (Section S7).
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

# maximum contextual confidence, directional form with O_l = 0 (Fig. S3)
kappa <- function(eps_Y, eps_l, f) (1 - eps_Y) * (1 - f) / ((1 - f) + eps_l)


# =============================================================================
# 5. CHECKING THE FUNCTION AGAINST THE SI
# =============================================================================
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

cat(sprintf("\n  kappa for THIS dataset = %.3f\n", kappa(EPS_Y, EPS_L, F_POS)))
cat("  No link can exceed this, however many replicates are added.\n")


# =============================================================================
# 6. THE RUNS
# =============================================================================
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
# evidence it is about to be combined with. R is nominal (5), matching the
# decision that every site is a replicate; co-occurrence stays in rho.

R_NOMINAL <- 5

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
         pi_r_deg = 1 - (1 - pi_l_deg)^R_NOMINAL) %>%
  select(focal_site, pollinator, plant, d_pol, d_pla, pi_l_deg, pi_r_deg)

ev <- ev %>% left_join(regional, by = c("focal_site", "pollinator", "plant"))

# ---- run scenarios ----

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
#   The practical difference is the 147 links that are predicted (Y = 1) yet
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
#              generalists start higher than two specialists, and the mean
#              pi_l of 0.09 is conservative against the observed 0.24.
#
# WHAT THE COMPARISONS ANSWER.
#   A-uniform            the baseline, and the direct counterpart of the
#                        deterministic categorisation
#   A-uniform vs A-degree   how much does ecological prior information move
#                        the answer, holding the model's role fixed?
#   A-uniform vs B-uniform  how much does keeping the model's score instead of
#                        its threshold move the answer, holding the prior
#                        fixed?
#   B-degree             both together, the most informed run

run_settings <- list(
  "A-uniform" = list(model_as = "evidence", pi_Y = 0.5,
                     pi_l = 0.5,         pi_r = 0.5),
  "A-degree"  = list(model_as = "evidence", pi_Y = 0.5,
                     pi_l = ev$pi_l_deg, pi_r = ev$pi_r_deg),
  "B-uniform" = list(model_as = "prior",    pi_Y = ev$q,
                     pi_l = 0.5,         pi_r = 0.5),
  "B-degree"  = list(model_as = "prior",    pi_Y = ev$q,
                     pi_l = ev$pi_l_deg, pi_r = ev$pi_r_deg)
)

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
    p1bar = P1,
    nu    = NU,                  # finite: heterogeneous p1, Section S7
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

# =============================================================================
# 7. THE MAP
# =============================================================================
# The richest site, one category, shaded by posterior. Same species, same
# ordering and same look as map_richest_site in the deterministic script, so
# the two can be laid side by side. An asterisk marks links the cameras found.

MAP_RUN      <- "B-uniform"          # the direct counterpart of the hard map is A-uniform. Switch to A-degree, B-uniform, B-degree if desirable
MAP_CATEGORY <- "possibly missing"

richest <- ev %>%
  filter(O_l == 1) %>%
  count(focal_site, sort = TRUE) %>%
  slice_head(n = 1) %>%
  pull(focal_site)

# order species by overall degree, as in the deterministic map, then keep only
# those present at this site
pol_levels <- links %>% count(higher_level, sort = TRUE) %>% pull(higher_level)
pla_levels <- links %>% count(lower_level,  sort = TRUE) %>% pull(lower_level)

map_df <- results %>%
  filter(run == MAP_RUN, focal_site == richest) %>%
  mutate(conf       = .data[[MAP_CATEGORY]],
         pollinator = factor(pollinator, levels = pol_levels),
         plant      = factor(plant,      levels = rev(pla_levels))) %>%
  droplevels()

# The colour scale stops at the highest posterior actually reached, not at 1.
# kappa caps every posterior, so a 0-to-1 scale would wash the figure out and
# quietly overstate confidence. The legend has to say where the top is.
CONF_MAX <- max(map_df$conf)

ital <- function(x) lapply(strsplit(x, " "),
                           function(y) bquote(italic(.(paste(y, collapse = " ")))))

map_bayes <- ggplot(map_df, aes(pollinator, plant, fill = conf)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(data = filter(map_df, camera == 1), aes(label = "*"),
            colour = "grey20", size = 5, vjust = 0.72) +
  scale_fill_gradient(low = "white", high = "#DBA040",
                      limits = c(0, CONF_MAX),
                      name = sprintf("P(%s)", MAP_CATEGORY)) +
  scale_x_discrete(labels = ital) +
  scale_y_discrete(labels = ital) +
  labs(x = "Pollinator", y = "Plant",
       title = sprintf("%s: confidence that each link is %s",
                       richest, MAP_CATEGORY),
       subtitle = sprintf("Scale stops at %.2f, the highest reached. * marks a camera record.",
                          CONF_MAX)) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 11, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y  = element_text(size = 11),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text  = element_text(size = 11),
    panel.grid   = element_blank()
  )

ggsave(file.path(OUT_DIR, "map_bayesian_richest_site.pdf"), map_bayes,
       width = 14, height = 7)

cat(sprintf("\nDone. Map for %s written; scale tops out at %.3f.\n",
            richest, CONF_MAX))
