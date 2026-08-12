## ---------------------------------------------------------------------------
## Main-text figure (fig:posterior_map): what the Bayesian reading adds to the
## deterministic taxonomy, without any appeal to accumulating replicates.
##   (a) posterior over categories for each evidence pattern, uniform prior
##   (b) the same under an informative prior; the maximum moves off the
##       deterministic cell for two of the four patterns
##   (c) how much prior scepticism is needed to move it, as a function of R
## Directional rates of the SI (Eq. erep) with R = 5 replicates.
## Only the four categories in which the model is right are shown; with
## pi_Y = 1/2 they hold exactly 1 - eps_Y = 80% of the posterior throughout.
## ---------------------------------------------------------------------------

library(ggplot2)
library(patchwork)

eY <- 0.20; el <- 0.30; f <- 0.05; rho <- 0.15
R0 <- 5                                   # replicates in (a) and (b)
pi_l <- 0.10; pi_r <- 0.20                # informative prior; pi_Y stays 1/2

p1    <- rho * (1 - el)
epsP  <- function(R) ((1 - p1)^R - (1 - rho)^R) / (1 - (1 - rho)^R)
epsM  <- function(R) 1 - (1 - f)^R

sig <- rbind(c(1,1,1), c(1,1,0), c(1,0,1), c(1,0,0),
             c(0,1,1), c(0,1,0), c(0,0,1), c(0,0,0))
nm  <- c("recurrent", "locally\nunique", "possibly\nmissing", "phantom",
         "model-\nelusive", "weakly\nsupported", "locally\nabsent",
         "possibly\nforbidden")

lik <- function(E, R) sapply(1:8, function(i) {
  fy <- if (E[1] == sig[i,1]) 1 - eY else eY
  fl <- if (sig[i,2] == 1) (if (E[2] == 1) 1 - el else el)
        else               (if (E[2] == 0) 1 - f  else f)
  fr <- if (sig[i,3] == 1) (if (E[3] == 1) 1 - epsP(R) else epsP(R))
        else               (if (E[3] == 1) epsM(R)     else 1 - epsM(R))
  fy * fl * fr
})
prior <- function(pY, pl, pr) {
  v <- c(pY, pl, pr)
  sapply(1:8, function(i) prod(v^sig[i,] * (1 - v)^(1 - sig[i,])))
}
post <- function(E, pri, R = R0) { L <- lik(E, R) * pri; L / sum(L) }

## ---- the 4 x 4 grids --------------------------------------------------------
rows <- 1:4                                  # evidence patterns with Y = 1
cols <- 1:4                                  # categories with z_Y = 1
ev_lab <- apply(sig[rows, ], 1, function(z) sprintf("(%d, %d, %d)", z[1], z[2], z[3]))

grid_for <- function(pri) {
  do.call(rbind, lapply(rows, function(i) {
    p <- post(sig[i, ], pri)
    data.frame(ev = ev_lab[i], cat = nm[cols], p = p[cols],
               deterministic = cols == i,
               map = cols == which.max(p), stringsAsFactors = FALSE)
  }))
}
g_uni <- grid_for(prior(0.5, 0.5, 0.5))
g_inf <- grid_for(prior(0.5, pi_l, pi_r))

lvl_ev  <- rev(ev_lab)                       # first pattern at the top
lvl_cat <- nm[cols]
fix <- function(d) { d$ev <- factor(d$ev, levels = lvl_ev)
                     d$cat <- factor(d$cat, levels = lvl_cat); d }
g_uni <- fix(g_uni); g_inf <- fix(g_inf)

heat <- function(d, ttl) {
  ggplot(d, aes(cat, ev)) +
    geom_tile(aes(fill = p), colour = "white", linewidth = 0.8) +
    geom_tile(data = subset(d, map), fill = NA, colour = "#111111",
              linewidth = 1.1, width = 0.86, height = 0.86) +
    geom_point(data = subset(d, deterministic), shape = 21, size = 1.9,
               stroke = 0.6, fill = "grey35", colour = "white",
               position = position_nudge(x = -0.34, y = 0.34)) +
    geom_text(aes(label = sprintf("%.2f", p),
                  colour = p > 0.50), size = 3, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20")) +
    scale_fill_gradient(low = "#EAF1FA", high = "#1B4F8C",
                        limits = c(0, 0.75), name = "posterior") +
    labs(title = ttl, x = NULL, y = expression(bold("evidence ") *
         bold("(") * bolditalic(Y) * bold(", ") * bolditalic(O)[bold(l)] *
         bold(", ") * bolditalic(O)[bold(r)] * bold(")"))) +
    coord_fixed() +
    theme_classic(base_size = 10) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          plot.title = element_text(size = 10, hjust = 0),
          legend.position = "none")
}

pa <- heat(g_uni, "(a) uniform prior")
pb <- heat(g_inf, bquote("(b) informative prior: " * pi[l] ~ "=" ~ .(pi_l) *
                         "," ~ pi[r] ~ "=" ~ .(pi_r)))

## ---- (c) the prior needed to overturn the label -----------------------------
Rs  <- 1:20
thr <- 1 / (1 + (1 - epsP(Rs)) / epsM(Rs))
dc  <- data.frame(R = Rs, thr = thr)

pc <- ggplot(dc, aes(R, thr)) +
  geom_ribbon(aes(ymin = 0, ymax = thr), fill = "#E07B39", alpha = 0.20) +
  geom_line(linewidth = 0.8, colour = "#E07B39") +
  annotate("point", x = R0, y = thr[R0], size = 2.2, colour = "#111111") +
  annotate("text", x = R0 + 0.6, y = thr[R0], label = "case in (b)",
           hjust = 0, vjust = 1.6, size = 3, colour = "grey20") +
  annotate("text", x = 1, y = 0.44, hjust = 0, size = 3, colour = "grey25",
           label = "the taxonomy's label\nsurvives the prior") +
  annotate("text", x = 19, y = 0.06, hjust = 1, size = 3, colour = "#9C4E13",
           label = "the prior overturns it") +
  scale_y_continuous(limits = c(0, 0.55),
                     labels = function(x) sprintf("%.1f", x)) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  labs(title = "(c) prior needed to overturn the label",
       x = expression("Replicates, " * italic(R)),
       y = expression("prior on regional presence, " * pi[r])) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0),
        axis.line = element_line(colour = "grey40", linewidth = 0.3),
        axis.ticks = element_line(colour = "grey40", linewidth = 0.3))

fig <- pa + pb + pc + plot_layout(widths = c(1, 1, 1.05))

## ---- reported values --------------------------------------------------------
cat(sprintf("eps_r+(%d) = %.3f   eps_r-(%d) = %.3f\n", R0, epsP(R0), R0, epsM(R0)))
for (nmv in c("uniform", "informative")) {
  d <- if (nmv == "uniform") g_uni else g_inf
  cat("\n", nmv, "\n")
  for (e in ev_lab) {
    s <- subset(d, ev == e)
    cat(sprintf("  %s  %s | MAP %s %.3f | deterministic %.3f\n", e,
        paste(sprintf("%.3f", s$p), collapse = " "),
        gsub("\n", " ", as.character(s$cat[s$map])), s$p[s$map], s$p[s$deterministic]))
  }
}
cat("\nthreshold pi_r:", sprintf("R=%d %.3f", c(1,3,5,10,20), thr[c(1,3,5,10,20)]), "\n")

ggsave("posterior_map.pdf", fig, width = 11.0, height = 3.6)
ggsave("posterior_map.png", fig, width = 11.0, height = 3.6, dpi = 300)
