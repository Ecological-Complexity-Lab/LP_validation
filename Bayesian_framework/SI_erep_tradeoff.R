# =============================================================================
# Supplementary Information - the replicate axis as an error-prone evidence
# source. Two figures.
#
# Model (same generative model as the replicate term of accumulation_curves.R):
#   a truly present link is realised in a replicate with probability rho and
#   then detected with the local sampling success (1 - eps_L), so
#       p1 = rho * (1 - eps_L),
#   and a truly absent link is spuriously recorded with probability
#       p0 = f.
#   Summarising R replicates as a single "seen at least once" bit gives two
#   directional error rates,
#       eps_rep_present(R) = (1 - p1)^R        false negative, falls with R
#       eps_rep_absent(R)  = 1 - (1 - p0)^R    false positive, rises with R
#   R is treated as continuous throughout.
#
#   erep_tradeoff.{pdf,png}     the tension between the two directional rates
#   evidential_weight.{pdf,png} the likelihood ratios carried by a detection
#                               and by an empty replicate set
# =============================================================================

library(ggplot2)
library(patchwork)

## ---- parameters -------------------------------------------------------------
eps_Y <- 0.20    # model error (epsilon_Y); taxonomy only, not these figures
eps_L <- 0.30    # local-miss error  (epsilon_local)
rho   <- 0.15    # realisation rate
f     <- 0.05    # false-detection rate

Rmax_tension <- 30              # x range of figure 1, panel (a)
Rmax_weight  <- 50              # x range of figure 2
rho_range    <- c(0.02, 0.90)   # x range of figure 1, panel (b)
rho_marks    <- c(0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.80)  # labelled points
## hand-placed x offsets for those labels; the wider one clears the rho = 0.15 line
rho_mark_dx  <- c(0.020, 0.055, 0.025, 0.020, 0.020, 0.020, 0.020)
rho_mark_dy  <- 0.035
f_levels     <- c(0.01, 0.02, 0.05, 0.10)   # one curve each, figure 2

fig_w <- 9.5     # inches
fig_h <- 4.6
fig_dpi <- 300

## ---- palette (ColorBrewer Dark2) --------------------------------------------
dark2 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
           "#66A61E", "#E6AB02", "#A6761D", "#666666")

## ---- the model --------------------------------------------------------------
p1_of <- function(rho, eps_L) rho * (1 - eps_L)   # per-replicate true detection
p0_of <- function(f) f                            # per-replicate false detection

## false negative after R replicates: a realisable link (L_r = 1) goes unrecorded
## only if it is unobserved in every replicate. Realisation is not conditioned on,
## because the truth being recovered is realisability, not realisation, so the rate
## carries the ecological failure to be realised as well as the failure to detect.
eps_present <- function(R, p1, rho = NULL) (1 - p1)^R
eps_absent  <- function(R, p0) 1 - (1 - p0)^R     # false positive after R replicates

## R at which the two directional rates are equal. eps_present falls
## monotonically from 1 and eps_absent rises monotonically from 0, so the
## difference is monotone and the root is unique.
crossing_R <- function(p1, p0, rho = NULL) {
  g <- function(R) eps_present(R, p1, rho) - eps_absent(R, p0)
  uniroot(g, interval = c(1, 1e4), tol = .Machine$double.eps^0.5)$root
}

## the common rate at that crossing = the lowest value the pointwise maximum
## of the two rates can take
best_rate <- function(p1, p0, rho = NULL) eps_present(crossing_R(p1, p0), p1)

p1 <- p1_of(rho, eps_L)
p0 <- p0_of(f)

cat(sprintf("Per-replicate detection probabilities\n"))
cat(sprintf("  p1 = rho * (1 - eps_L) = %.3f * %.3f = %.4f\n", rho, 1 - eps_L, p1))
cat(sprintf("  p0 = f                 = %.4f\n\n", p0))

## ---- numerical checks -------------------------------------------------------
## Values quoted in the SI text. tol = 5e-4 means agreement to three decimal
## places; the crossing values of R are quoted to two, hence tol = 5e-3.
chk <- function(quantity, computed, expected, tol = 5e-4) {
  data.frame(quantity = quantity, expected = expected, computed = computed,
             tol = tol, ok = abs(computed - expected) < tol,
             stringsAsFactors = FALSE)
}

R_star    <- crossing_R(p1, p0, rho)
rate_star <- eps_present(R_star, p1, rho)
mx        <- optimize(function(R) pmax(eps_present(R, p1, rho), eps_absent(R, p0)),
                      interval = c(1, Rmax_tension), tol = 1e-8)

panel_b_expected <- data.frame(
  rho  = c(0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.80),
  R    = c(16.13, 11.31, 8.98, 7.54, 5.77, 3.94, 2.55),
  rate = c(0.563, 0.440, 0.369, 0.321, 0.256, 0.183, 0.123)
)

w_det   <- function(R, p1, p0) (1 - (1 - p1)^R) / (1 - (1 - p0)^R)
w_empty <- function(R, p1, p0) ((1 - p1) / (1 - p0))^R
R_probe <- c(1, 5, 10, 20, 50)

checks <- do.call(rbind, c(
  lapply(seq_along(c(1, 5, 10)), function(i) {
    R <- c(1, 5, 10)[i]
    chk(sprintf("eps_rep^present, R = %d", R), eps_present(R, p1, rho),
        c(0.895, 0.574, 0.330)[i])
  }),
  lapply(seq_along(c(1, 5, 10)), function(i) {
    R <- c(1, 5, 10)[i]
    chk(sprintf("eps_rep^absent, R = %d", R), eps_absent(R, p0),
        c(0.050, 0.226, 0.401)[i])
  }),
  list(
    chk("crossing, R",                R_star,     8.98,  tol = 5e-3),
    chk("crossing, common rate",      rate_star,  0.369),
    chk("min of pointwise max, R",    mx$minimum, 8.98,  tol = 5e-3),
    chk("min of pointwise max, rate", mx$objective, 0.369)
  ),
  lapply(seq_len(nrow(panel_b_expected)), function(i) {
    r <- panel_b_expected$rho[i]
    chk(sprintf("panel (b) rho = %.2f, R", r),
        crossing_R(p1_of(r, eps_L), p0, r), panel_b_expected$R[i], tol = 5e-3)
  }),
  lapply(seq_len(nrow(panel_b_expected)), function(i) {
    r <- panel_b_expected$rho[i]
    chk(sprintf("panel (b) rho = %.2f, rate", r),
        best_rate(p1_of(r, eps_L), p0, r), panel_b_expected$rate[i])
  }),
  lapply(seq_along(R_probe), function(i) {
    chk(sprintf("detection weight, f = 0.05, R = %d", R_probe[i]),
        w_det(R_probe[i], p1, 0.05),
        c(2.100, 1.882, 1.670, 1.389, 1.079)[i])
  }),
  lapply(seq_along(R_probe), function(i) {
    chk(sprintf("empty-set weight, f = 0.05, R = %d", R_probe[i]),
        w_empty(R_probe[i], p1, 0.05),
        c(0.942, 0.742, 0.551, 0.303, 0.051)[i])
  })
))

checks_print <- checks
checks_print$computed <- round(checks_print$computed, 4)
cat("Numerical checks\n")
print(checks_print, row.names = FALSE)
cat("\n")

if (any(!checks$ok)) {
  bad <- checks[!checks$ok, ]
  stop(sprintf(
    "%d numerical check(s) failed:\n%s",
    nrow(bad),
    paste(sprintf("  %-38s expected %.3f, computed %.6f (tol %.0e)",
                  bad$quantity, bad$expected, bad$computed, bad$tol),
          collapse = "\n")), call. = FALSE)
}

# =============================================================================
# FIGURE 1 - erep_tension
# =============================================================================

## ---- panel (a): the two directional rates and their maximum -----------------
R_grid <- seq(1, Rmax_tension, length.out = 601)

dat_a <- data.frame(
  R      = rep(R_grid, 3),
  series = factor(rep(c("present", "absent", "max"), each = length(R_grid)),
                  levels = c("present", "absent", "max")),
  value  = c(eps_present(R_grid, p1, rho),
             eps_absent(R_grid, p0),
             pmax(eps_present(R_grid, p1, rho), eps_absent(R_grid, p0)))
)

pal_a <- c(present = dark2[1], absent = dark2[2], max = "grey55")
lty_a <- c(present = "solid",  absent = "solid",  max = "dashed")
lab_a <- expression(epsilon[r]^{"+"} ~ "(false negative)",
                    epsilon[r]^{"-"} ~ "(false positive)",
                    "maximum of the two")

## The maximum coincides with whichever directional rate is currently larger, so
## it is drawn first and wider: it reads as a grey dashed envelope flanking the
## two solid curves rather than being hidden underneath them. override.aes below
## puts every legend key back at a common line width.
p_a <- ggplot() +
  geom_line(data = subset(dat_a, series == "max"),
            aes(R, value, colour = series, linetype = series), linewidth = 2.1) +
  geom_line(data = subset(dat_a, series != "max"),
            aes(R, value, colour = series, linetype = series), linewidth = 0.9) +
  ## the crossing, with a leader up to the annotation
  annotate("segment", x = R_star, xend = R_star,
           y = rate_star + 0.03, yend = 0.875,
           linetype = "dotted", colour = "grey40", linewidth = 0.35) +
  annotate("point", x = R_star, y = rate_star, size = 2.6, colour = "grey20") +
  annotate("text", x = R_star, y = 0.925,
           label = sprintf("R = %.2f,  rate = %.3f", R_star, rate_star),
           size = 3.1, colour = "grey20") +
  scale_colour_manual(name = NULL, values = pal_a,
                      breaks = names(pal_a), labels = lab_a) +
  scale_linetype_manual(name = NULL, values = lty_a,
                        breaks = names(lty_a), labels = lab_a) +
  guides(colour = guide_legend(override.aes = list(linewidth = 0.9))) +
  scale_x_continuous(breaks = seq(0, Rmax_tension, 5),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = expression("Number of replicates, " * italic(R)),
       y = "Error rate") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.text      = element_text(size = 10))

## ---- panel (b): the best the replicate axis can do, as a function of rho ----
rho_grid <- seq(rho_range[1], rho_range[2], length.out = 401)

dat_b <- data.frame(
  rho  = rho_grid,
  R    = vapply(rho_grid, function(r) crossing_R(p1_of(r, eps_L), p0, r), numeric(1))
)
dat_b$rate <- mapply(function(R, r) eps_present(R, p1_of(r, eps_L), r), dat_b$R, dat_b$rho)

marks_b <- data.frame(
  rho = rho_marks,
  R   = vapply(rho_marks, function(r) crossing_R(p1_of(r, eps_L), p0, r), numeric(1))
)
marks_b$rate  <- mapply(function(R, r) eps_present(R, p1_of(r, eps_L), r), marks_b$R, marks_b$rho)
marks_b$label <- sprintf("R = %.1f", marks_b$R)
marks_b$lab_x <- marks_b$rho  + rho_mark_dx
marks_b$lab_y <- marks_b$rate + rho_mark_dy

p_b <- ggplot(dat_b, aes(rho, rate)) +
  geom_vline(xintercept = rho, linetype = "dashed", colour = "grey55") +
  geom_line(colour = dark2[3], linewidth = 1) +
  geom_point(data = marks_b, size = 2, colour = "grey20") +
  geom_text(data = marks_b, aes(lab_x, lab_y, label = label), hjust = 0,
            size = 3, colour = "grey20") +
  scale_x_continuous(breaks = seq(0, 0.9, 0.2),
                     expand = expansion(mult = c(0.02, 0.10))) +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 0.8, 0.2),
                     expand = expansion(mult = c(0.01, 0.08))) +
  labs(x = expression(rho), y = "Best achievable error rate") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

fig1 <- (p_a | p_b) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom",
        plot.tag = element_text(size = 12, face = "bold"))

print(fig1)
ggsave("SI_erep_tradeoff.pdf", fig1, width = fig_w, height = fig_h)
ggsave("SI_erep_tradeoff.png", fig1, width = fig_w, height = fig_h, dpi = fig_dpi)

# =============================================================================
# FIGURE 2 - evidential_weight
# =============================================================================
R_grid2 <- seq(1, Rmax_weight, length.out = 501)

dat_w <- expand.grid(R = R_grid2, f = f_levels)
dat_w$detection <- w_det(dat_w$R, p1, dat_w$f)
dat_w$empty     <- w_empty(dat_w$R, p1, dat_w$f)
dat_w$f_lab     <- factor(sprintf("%.2f", dat_w$f),
                          levels = sprintf("%.2f", f_levels))

pal_f <- setNames(dark2[seq_along(f_levels)], levels(dat_w$f_lab))
lab_f <- expression("False-detection rate, " * italic(f))

weight_theme <- list(
  scale_colour_manual(name = lab_f, values = pal_f),
  scale_x_continuous(breaks = c(1, seq(10, Rmax_weight, 10)),
                     expand = expansion(mult = c(0.02, 0.02))),
  labs(x = expression("Number of replicates, " * italic(R))),
  theme_bw(base_size = 12),
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom")
)

p_w1 <- ggplot(dat_w, aes(R, detection, colour = f_lab)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 1) +
  labs(y = "Weight of a detection") +
  weight_theme

p_w2 <- ggplot(dat_w, aes(R, empty, colour = f_lab)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 1) +
  scale_y_log10() +
  labs(y = "Weight of an empty replicate set") +
  weight_theme

fig2 <- (p_w1 | p_w2) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom",
        plot.tag = element_text(size = 12, face = "bold"))

print(fig2)
ggsave("SI_evidential_weight.pdf", fig2, width = fig_w, height = fig_h)
ggsave("SI_evidential_weight.png", fig2, width = fig_w, height = fig_h, dpi = fig_dpi)

cat("Wrote erep_tradeoff.{pdf,png} and evidential_weight.{pdf,png}\n")
