## ---------------------------------------------------------------------------
## SI figure (fig:model_score): the model score as prior against the model as
## thresholded evidence. Evidence on the other two axes is fixed at O_l = 0 and
## detection in every one of R = 5 replicates.
##   pi_Y = q      : the score enters the prior, no model factor in the likelihood
##   thresholded   : Y = 1[q > tau] enters as evidence with rate eps_Y, pi_Y = 1/2
## ---------------------------------------------------------------------------

library(ggplot2)

## Figures are written to Bayesian_framework/bayesian_figures. Resolved here so
## the script works whether it is run from the project root or from inside
## Bayesian_framework/, and the folder is created if it is missing.
fig_dir <- if (dir.exists("Bayesian_framework"))
             file.path("Bayesian_framework", "bayesian_figures") else "bayesian_figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

eY <- 0.20; el <- 0.30; f <- 0.05; rho <- 0.15
R <- 5; tau <- 0.5

p1 <- rho * (1 - el); eP <- (1 - p1)^R; eM <- 1 - (1 - f)^R
sig <- rbind(c(1,1,1), c(1,1,0), c(1,0,1), c(1,0,0),
             c(0,1,1), c(0,1,0), c(0,0,1), c(0,0,0))
nm  <- c("recurrent","locally unique","possibly missing","phantom",
         "model-elusive","weakly supported","locally absent","possibly forbidden")

## local and replicate factors for O_l = 0, O_r = 1 (identical in the two halves)
lr <- sapply(1:8, function(i) {
  fl <- if (sig[i,2] == 1) el else 1 - f
  fr <- if (sig[i,3] == 1) 1 - eP else eM
  fl * fr
})

## model as prior: pi_Y = q, no model factor in the likelihood
as_prior <- function(q) {
  pri <- sapply(1:8, function(i) ifelse(sig[i,1] == 1, q, 1 - q) * 0.25)
  L <- lr * pri; L / sum(L)
}
## model as evidence: Y = 1[q > tau], pi_Y = 1/2, model factor eps_Y
as_evidence <- function(q) {
  Y <- as.integer(q > tau)
  L <- lr * sapply(1:8, function(i) if (sig[i,1] == Y) 1 - eY else eY) * 0.125
  L / sum(L)
}

qs <- seq(0.001, 0.999, length.out = 999)
shown <- 3                                        # possibly missing
dat <- do.call(rbind, lapply(shown, function(i) rbind(
  data.frame(q = qs, p = sapply(qs, function(q) as_prior(q)[i]),
             cat = nm[i], how = "score as prior"),
  data.frame(q = qs, p = sapply(qs, function(q) as_evidence(q)[i]),
             cat = nm[i], how = "thresholded"))))
dat$cat <- factor(dat$cat, levels = nm[shown])
dat$how <- factor(dat$how, levels = c("score as prior", "thresholded"))
dat$grp <- interaction(dat$cat, dat$how)
## break the step so the vertical jump is not drawn
dat$grp <- interaction(dat$grp, dat$q > tau)

pal <- setNames("#E07B39", nm[shown])
## the vertical jump of the step, drawn so it reads as a step function
jump <- data.frame(q = c(tau, tau),
                   p = c(as_evidence(tau - 1e-6)[3], as_evidence(tau + 1e-6)[3]))
meet <- data.frame(q = c(eY, 1 - eY),
                   p = c(as_prior(eY)[3], as_prior(1 - eY)[3]))

p <- ggplot(dat, aes(q, p, colour = cat, linetype = how, group = grp)) +
  geom_vline(xintercept = tau, colour = "grey75", linewidth = 0.35) +
  geom_line(linewidth = 0.8) +
  geom_line(data = jump, aes(q, p), inherit.aes = FALSE, linetype = "22",
            colour = pal[[1]], linewidth = 0.8) +
  geom_point(data = meet, aes(q, p), inherit.aes = FALSE, size = 2,
             colour = "grey20") +
  annotate("text", x = eY, y = as_prior(eY)[3], label = "q == epsilon[Y]",
           hjust = 1.15, vjust = -0.4, size = 3, colour = "grey25", parse = TRUE) +
  annotate("text", x = 1 - eY, y = as_prior(1 - eY)[3],
           label = "q == 1 - epsilon[Y]",
           hjust = 1.1, vjust = -0.6, size = 3, colour = "grey25", parse = TRUE) +
  annotate("text", x = tau, y = 0.05, label = "threshold", angle = 90,
           vjust = -0.5, hjust = 0, size = 3, colour = "grey55") +
  scale_colour_manual(values = pal, guide = "none") +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  scale_y_continuous(limits = c(0, 0.55),
                     labels = function(x) paste0(round(100 * x), "%")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = expression("Calibrated model score, " * italic(q)),
       y = expression(P(italic("possibly missing") ~ "|" ~ E))) +
  theme_classic(base_size = 10) +
  theme(axis.line = element_line(colour = "grey40", linewidth = 0.3),
        axis.ticks = element_line(colour = "grey40", linewidth = 0.3),
        legend.position = "bottom", legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"))

cat(sprintf("slope of the ramp for possibly missing: %.4f per unit q\n", as_prior(1)[3]))
for (q in c(0.05, 0.20, 0.49, 0.51, 0.80, 0.95))
  cat(sprintf("  q=%.2f  prior %.3f / %.3f   thresholded %.3f / %.3f   (pm / loc.absent)\n",
      q, as_prior(q)[3], as_prior(q)[7], as_evidence(q)[3], as_evidence(q)[7]))

ggsave(file.path(fig_dir, "SI_prior_q.pdf"), p, width = 5.4, height = 3.6)
ggsave(file.path(fig_dir, "SI_prior_q.png"), p, width = 5.4, height = 3.6, dpi = 300)
