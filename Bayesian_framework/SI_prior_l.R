## ---------------------------------------------------------------------------
## SI figure (fig:degree_prior): what a degree-based prior does to a fixed piece
## of evidence. The evidence is the one used throughout, E = (1, 0, 1): the link
## is predicted, not recorded locally, recorded in the replicates. Both priors
## rise with the pair's degree, pi_l with local degree and pi_r with regional
## degree (Section prior). The mapping from degree to pi is an illustration,
## not a derivation.
##   (a) pi_l = 0.9 pi_r, varying pi_r: how generalist the pair is
##   (b) pi_r = 0.9 fixed, varying c = pi_l/pi_r: how much of the regional
##       expectation transfers to the focal site
## Only the four categories in which the model is right are drawn; with
## pi_Y = 1/2 they hold exactly 1 - eps_Y = 80% of the posterior throughout.
## ---------------------------------------------------------------------------

library(ggplot2)
library(patchwork)

eY <- 0.20; el <- 0.30; f <- 0.05; rho <- 0.15
R  <- 5

p1 <- rho * (1 - el)
eP <- (1 - p1)^R          # eps_r^+(R), Eq. erep
eM <- 1 - (1 - f)^R       # eps_r^-(R)

sig <- rbind(c(1,1,1), c(1,1,0), c(1,0,1), c(1,0,0),
             c(0,1,1), c(0,1,0), c(0,0,1), c(0,0,0))
nm  <- c("recurrent", "locally unique", "possibly missing", "phantom",
         "model-elusive", "weakly supported", "locally absent",
         "possibly forbidden")

## posterior for E = (1, 0, 1) under Eq. prior with pi_Y = 1/2
post <- function(pl, pr) {
  pri <- sapply(1:8, function(i)
    prod(c(0.5, pl, pr)^sig[i, ] * (1 - c(0.5, pl, pr))^(1 - sig[i, ])))
  L <- sapply(1:8, function(i) {
    fy <- if (sig[i,1] == 1) 1 - eY else eY
    fl <- if (sig[i,2] == 1) el else 1 - f      # O_l = 0, directional
    fr <- if (sig[i,3] == 1) 1 - eP else eM     # O_r = 1
    fy * fl * fr
  }) * pri
  L / sum(L)
}

shown <- 1:4                                     # the z_Y = 1 block
pal <- setNames(c("#2A78D6", "#2A9D8F", "#E07B39", "#7C5CBF"), nm[shown])

grid <- seq(0.02, 0.95, length.out = 400)

build <- function(x, pl_of, pr_of, ttl, xlab) {
  P <- sapply(x, function(v) post(pl_of(v), pr_of(v)))
  d <- do.call(rbind, lapply(shown, function(i)
    data.frame(pr = x, p = P[i, ], cat = nm[i])))
  d$cat <- factor(d$cat, levels = nm[shown])
  ## where the leading category changes
  lead <- nm[apply(P[shown, ], 2, which.max)]
  brk  <- x[which(diff(as.integer(factor(lead, levels = nm[shown]))) != 0)]
  segs <- c(min(x), brk, max(x))
  labs <- data.frame(
    x = (head(segs, -1) + tail(segs, -1)) / 2, y = 0.94,
    lab = lead[sapply((head(segs, -1) + tail(segs, -1)) / 2,
                      function(v) which.min(abs(x - v)))])
  ggplot(d, aes(pr, p, colour = cat)) +
    { if (length(brk)) geom_vline(xintercept = brk, linetype = "dotted",
                                  colour = "grey60", linewidth = 0.35) } +
    geom_line(linewidth = 0.8) +
    geom_text(data = labs, aes(x, y, label = lab), inherit.aes = FALSE,
              size = 2.9, colour = "grey30", fontface = "italic") +
    scale_colour_manual(values = pal, name = NULL) +
    scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(round(100*x), "%")) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(title = ttl, x = xlab, y = expression(P(C ~ "|" ~ E))) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(size = 10, hjust = 0),
          axis.line = element_line(colour = "grey40", linewidth = 0.3),
          axis.ticks = element_line(colour = "grey40", linewidth = 0.3),
          legend.position = "bottom")
}

cA <- 0.9; prB <- 0.9
pa <- build(grid, function(v) cA * v, function(v) v,
            bquote("(a) how generalist the pair is (" * pi[l] == .(cA) * pi[r] * ")"),
            expression("Prior that the link is realisable regionally, " * pi[r]))
pb <- build(seq(0.02, 1, length.out = 400),
            function(v) v * prB, function(v) prB,
            bquote("(b) how much transfers to the site (" * pi[r] == .(prB) * ")"),
            expression(atop("Local relative to regional, " * pi[l]/pi[r],
                            scriptstyle("rare, non-overlapping" %->%
                                        "abundant, overlapping")))) +
  labs(y = NULL)

fig <- pa + pb + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

## ---- reported values --------------------------------------------------------
for (k in c(cA)) {
  cat(sprintf("\npanel (a): pi_l = %.2f pi_r\n", k))
  cat(sprintf("%6s %11s %11s %11s %11s  leader\n", "pi_r", nm[1], nm[2], nm[3], nm[4]))
  for (d in c(0.05, 0.10, 0.30, 0.50, 0.70, 0.90)) {
    p <- post(k * d, d)
    cat(sprintf("%6.2f %11.3f %11.3f %11.3f %11.3f  %s\n",
                d, p[1], p[2], p[3], p[4], nm[which.max(p[1:4])]))
  }
  P <- sapply(grid, function(d) post(k * d, d))
  lead <- apply(P[1:4, ], 2, which.max)
  b <- grid[which(diff(lead) != 0)]
  cat("  leader changes at pi_r =", sprintf("%.2f", b), "\n")
}
cat(sprintf("\npanel (b): pi_r = %.2f, varying c = pi_l/pi_r\n", prB))
cat(sprintf("%6s %11s %11s %11s %11s  leader\n", "c", nm[1], nm[2], nm[3], nm[4]))
for (cc in c(0.05, 0.25, 0.50, 0.75, 0.90, 1.00)) {
  p <- post(cc * prB, prB)
  cat(sprintf("%6.2f %11.3f %11.3f %11.3f %11.3f  %s\n", cc, p[1], p[2], p[3], p[4],
              nm[which.max(p[1:4])]))
}
cs <- seq(0.02, 1, length.out = 400)
lb <- apply(sapply(cs, function(v) post(v * prB, prB))[1:4, ], 2, which.max)
cat("  leader changes at c =", sprintf("%.2f", cs[which(diff(lb) != 0)]), "\n")

ggsave("SI_prior_l.pdf", fig, width = 7.6, height = 3.4)
ggsave("SI_prior_l.png", fig, width = 7.6, height = 3.4, dpi = 300)
