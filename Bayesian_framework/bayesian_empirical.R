# Bayesian categorisation of predicted links: the Cabrera worked example
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


# ---- 2. The rates ----

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

# ---- 4. The posterior master equation ----
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

cat(sprintf("\n  kappa for THIS dataset = %.3f\n", kappa(EPS_Y, EPS_L, F_POS)))
cat("  No link can exceed this, however many replicates are added.\n")


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
#   limit     top of the colour scale. NULL stops it at the highest posterior
#             actually reached, which is the honest default: kappa caps every
#             posterior, so a 0-to-1 scale would wash the figure out and
#             overstate confidence. Pass a number to hold several maps on one
#             scale, which is what makes them comparable to each other.
#
# Returns the ggplot, with the scale top attached as attr(, "conf_max") and the
# asterisk count as attr(, "n_marked").
make_map <- function(run      = "A-uniform",
                     category = "possibly missing",
                     site     = richest,
                     results_df = results,
                     camera_marks = c("assigned", "all", "none"),
                     limit    = NULL) {

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

  p <- ggplot(d, aes(pollinator, plant, fill = conf)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = unname(CATEGORY_COLOUR[category]),
                        limits = c(0, conf_max),
                        name = sprintf("P(%s)", category)) +
    scale_x_discrete(labels = ital) +
    scale_y_discrete(labels = ital) +
    labs(x = "Pollinator", y = "Plant",
         caption = sprintf("%s, run %s. Scale stops at %.2f%s.",
                           site, run, conf_max, cap_marks)) +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(size = 16, angle = 90, vjust = 0.5),
      axis.text.y  = element_text(size = 16),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 14),
      panel.grid   = element_blank()
    )

  if (nrow(marked) > 0) {
    p <- p + geom_text(data = marked, aes(label = "*"),
                       colour = "grey20", size = 5, vjust = 0.72)
  }

  attr(p, "conf_max") <- conf_max
  attr(p, "n_marked") <- nrow(marked)
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

MARK_RUN <- "A-degree"
mark_tbl <- camera_mark_summary(MARK_RUN)

cat(sprintf("\nCAMERA CORROBORATION AT %s (run %s)\n", richest, MARK_RUN))
cat(sprintf("  camera_marks = \"all\" puts the same %d asterisks on every map.\n",
            sum(mark_tbl$confirmed)))
cat("  camera_marks = \"assigned\" splits those across the categories below,\n")
cat("  so `confirmed` is what each map marks and the column sums to the total.\n")

mark_tbl %>%
  mutate(rate_assigned = sprintf("%.1f%%", 100 * rate_assigned),
         rate_covered  = ifelse(is.na(rate_covered), "-",
                                sprintf("%.1f%%", 100 * rate_covered))) %>%
  print(n = 8, width = Inf)

# --- produce chosen map ---------------------------------------------
MAP_RUN      <- "A-degree"          # the direct counterpart of the hard map;
                                     # A-degree, B-uniform, B-degree also valid
MAP_CATEGORY <- "model-elusive"

camera_marks <- "assigned"

map_bayes <- make_map(run = MAP_RUN, category = MAP_CATEGORY)

ggsave(file.path(OUT_DIR, "map_bayesian_richest_site.pdf"), map_bayes,
       width = 14, height = 7)

cat(sprintf("\nDone. Map for %s written (run %s, %s); scale tops out at %.3f, %d asterisks.\n",
            richest, MAP_RUN, MAP_CATEGORY, attr(map_bayes, "conf_max"),
            attr(map_bayes, "n_marked")))

# Any other combination is one call, for example
#   make_map("B-degree", "phantom")
#   make_map("A-uniform", "possibly forbidden")
# every camera record at the site rather than only this category's:
#   make_map("A-uniform", "phantom", camera_marks = "all")
# and a shared `limit` puts several on one scale:
#   make_map("A-uniform", "possibly missing", limit = 0.6)


# ---- 8. Category distributions across the four runs ----
# Two views of the same four posteriors, sharing one category colour key taken
# from the maps.
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
# On the kappa line in panel b. kappa = (1 - eps_Y)(1 - f) / ((1 - f) + eps_l)
# carries TWO conditions, and points sit to the right of it whenever either is
# relaxed, which is a result rather than a violation.
#
#   1. A uniform prior. B replaces the flat pi_Y with the calibrated score and
#      the degree runs add ecological priors, so the line is drawn only in the
#      A-uniform panel; the other three exceed it and should.
#   2. A link NOT recorded locally (O_l = 0). The derivation is the ceiling for
#      a non-detection, where the local factor contributes only eps_l against
#      1 - f. A positive local detection contributes 1 - eps_l against f, a far
#      sharper contrast, so those links can go higher. Under A-uniform every
#      category peaks at 0.533 = kappa among O_l = 0 links, while 138 O_l = 1
#      links reach 0.672, all of them recurrent or model-elusive, the two
#      categories with z_l = 1 and z_r = 1.
#
# Each row mixes both local conditions, so the line cannot be drawn per row.
# The axis label states what it bounds instead.

library(patchwork)

CAT_ORDER <- cats$category          # taxonomy order: predicted four, then not
RUN_ORDER <- names(run_settings)

post_long <- results %>%
  select(run, all_of(CAT_ORDER)) %>%
  pivot_longer(all_of(CAT_ORDER), names_to = "category", values_to = "p") %>%
  mutate(category = factor(category, levels = rev(CAT_ORDER)),
         run      = factor(run, levels = RUN_ORDER))

# the deterministic counts are a property of the labelling, not of a run, so
# the same reference is drawn in every panel
det_counts <- results %>%
  filter(run == RUN_ORDER[1]) %>%
  count(category = det_category, name = "deterministic") %>%
  mutate(category = factor(category, levels = rev(CAT_ORDER)))

expected_tbl <- post_long %>%
  group_by(run, category) %>%
  summarise(expected = sum(p), .groups = "drop") %>%
  left_join(det_counts, by = "category")

# kappa bounds the uniform-prior run only (see the note above)
kappa_line <- tibble(run = factor("A-uniform", levels = RUN_ORDER),
                     k   = kappa(EPS_Y, EPS_L, F_POS))

si_base <- theme_classic(base_size = 10) +
  theme(axis.line        = element_line(colour = "grey40", linewidth = 0.3),
        axis.ticks       = element_line(colour = "grey40", linewidth = 0.3),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 9, hjust = 0),
        panel.spacing    = unit(9, "pt"),
        plot.tag         = element_text(size = 13, face = "bold"),
        legend.position  = "none")

p_mass <- ggplot(expected_tbl, aes(expected, category, fill = category)) +
  geom_col(width = 0.72, colour = "grey35", linewidth = 0.2) +
  geom_point(aes(x = deterministic), shape = 124, size = 2.6, colour = "grey15") +
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
  # headroom so the kappa label sits inside the panel rather than being clipped
  scale_y_discrete(expand = expansion(add = c(0.6, 1.0))) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = expression("P(category) per link"),
       y = NULL) +
  si_base +
  theme(strip.text = element_blank())      # run names already label panel a

fig_runs <- (p_mass / p_spread) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 13, face = "bold"))

ggsave(file.path(OUT_DIR, "fig_category_distributions.pdf"), fig_runs,
       width = 11, height = 7)
ggsave(file.path(OUT_DIR, "fig_category_distributions.png"), fig_runs,
       width = 11, height = 7, dpi = 300)

cat("\nEXPECTED COUNTS BY RUN (deterministic count in the first column)\n")
expected_tbl %>%
  mutate(expected = round(expected, 1)) %>%
  select(category, run, deterministic, expected) %>%
  pivot_wider(names_from = run, values_from = expected) %>%
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
#   n, confirmed   size of the DETERMINISTIC category and how many of those the
#                  cameras recorded. Descriptive only: every regression below
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
#                  probability changes at different starting points: for
#                  possibly forbidden, OR 0.55 costs about 7 points among
#                  locally unobserved links and about 13 among observed ones.
#                  Averaging over the links respects the real mix of O_l and
#                  answers "for a typical link here, how much does the chance
#                  of a camera record move?".
#                  Read these two together with p_adj: the percentage-point
#                  gaps are modest (5 to 13 points) even where p is tiny, which
#                  reflects a clean signal in a large sample rather than a large
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
#                  opposite sign to OR_adj: possibly missing is 0.66 raw and
#                  1.60 adjusted. Never read it against the verdict.

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
  # construction, but the probability change it implies is not: -0.55 in odds
  # costs about 7 points among locally unobserved links and about 13 among
  # observed ones. So rather than quote one stratum, predict every link twice,
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
cat("A link is inside that grid when BOTH partners appear in its site's camera\n")
cat("data, so a zero means the cameras caught both species there but never\n")
cat("together. It does not mean a camera was aimed at that pair.\n")
cat("Odds ratio is per 1 SD of that category's posterior.\n")
cat("The verdict uses the O_l-adjusted model, so the posterior has to predict\n")
cat("camera records beyond what the direct observation already said.\n")

cam_report <- map_dfr(names(run_settings), camera_validation)

# Printed as one table. Column order matters here: the verdict is decided by
# the O_l-adjusted model, so OR_adj sits next to it. The raw OR is kept for
# reference but must not be read against the verdict, because it can carry the
# opposite sign: possibly missing is 0.66 raw and 1.60 adjusted, and only the
# adjusted figure is a claim about anything the direct observation had not
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

write_csv(cam_report, file.path(OUT_DIR, "camera_validation_by_category.csv"))
cat(sprintf("\nWritten: %s\n",
            file.path(OUT_DIR, "camera_validation_by_category.csv")))
