## ---------------------------------------------------------------------------
## SI figure: heterogeneous links. p_1 is drawn once per link from Beta(a,b)
## with mean pbar_1 = rho(1-eps_l) and concentration nu = a+b, and held fixed
## across that link's replicates, so the replicate count is beta-binomial.
## Evidence: Y = 1, O_l = 0. Directional local axis (Eq. likelihood_dir).
## ---------------------------------------------------------------------------

library(ggplot2)
library(patchwork)

eY   <- 0.20
el   <- 0.30
f    <- 0.05
rho  <- 0.15
pbar <- rho * (1 - el)          # 0.105
nus  <- c(100, 40, 10)          # concentrations shown, plus the common rate
Rs   <- 1:20

## the eight categories as (z_Y, z_l, z_r)
sig <- rbind(c(1,1,1), c(1,1,0), c(1,0,1), c(1,0,0),
             c(0,1,1), c(0,1,0), c(0,0,1), c(0,0,0))
rownames(sig) <- c("recurrent","locally unique","possibly missing","phantom",
                   "model-elusive","weakly-supported","locally absent",
                   "possibly forbidden")
PM <- 3; PH <- 4

## replicate factor for a category with z_r = 1, binomial coefficient dropped
rep1 <- function(n, R, nu) {
  if (is.na(nu)) pbar^n * (1 - pbar)^(R - n)          # common rate
  else exp(lbeta(pbar * nu + n, (1 - pbar) * nu + R - n) -
           lbeta(pbar * nu, (1 - pbar) * nu))          # beta-binomial
}

post <- function(R, n, nu, which) {
  L <- numeric(8)
  for (i in 1:8) {
    lY <- if (sig[i,1] == 1) 1 - eY else eY
    lL <- if (sig[i,2] == 1) el else 1 - f            # O_l = 0, directional
    pr <- if (sig[i,3] == 1) rep1(n, R, nu) else f^n * (1 - f)^(R - n)
    L[i] <- lY * lL * pr
  }
  (L / sum(L))[which]
}

kappa <- (1 - eY) * (1 - f) / ((1 - f) + el)          # 0.608

lev  <- c("common rate", paste0("nu == ", nus))
cols <- setNames(c("#4D4D4D", "#2A78D6", "#E07B39", "#7C5CBF"), lev)
labs <- setNames(c("common rate", expression(nu == 100), expression(nu == 40), expression(nu == 10)), lev)

## ---- (a) the distribution of p_1 across links ------------------------------
gridp <- seq(1e-4, 0.5, length.out = 800)
dens <- do.call(rbind, lapply(seq_along(nus), function(k) {
  data.frame(level = lev[k + 1], p = gridp,
             d = dbeta(gridp, pbar * nus[k], (1 - pbar) * nus[k]))
}))
dens$level <- factor(dens$level, levels = lev)

pa <- ggplot(dens, aes(p, d, colour = level)) +
  geom_vline(xintercept = pbar, colour = cols[1], linewidth = 0.7) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = cols, labels = labs, drop = FALSE) +
  labs(title = "(a) detection rate across links",
       x = expression(p[1]), y = "density")

## ---- (b) detected in every replicate, (c) empty in every replicate ---------
curves <- do.call(rbind, lapply(c(NA, nus), function(nu) {
  k <- if (is.na(nu)) lev[1] else paste0("nu == ", nu)
  rbind(data.frame(level = k, panel = "b", R = Rs,
                   p = sapply(Rs, function(R) post(R, R, nu, PM))),
        data.frame(level = k, panel = "c", R = Rs,
                   p = sapply(Rs, function(R) post(R, 0, nu, PM))))
}))
curves$level <- factor(curves$level, levels = lev)

pct <- function(x) paste0(round(100 * x), "%")
mk <- function(pn, ttl, yl, xmax)
  ggplot(subset(curves, panel == pn & R <= xmax), aes(R, p, colour = level)) +
    geom_hline(yintercept = kappa, linetype = "22", colour = "grey25",
               linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = cols, labels = labs, drop = FALSE) +
    scale_y_continuous(limits = c(0, 0.72), labels = pct) +
    scale_x_continuous(breaks = scales::pretty_breaks(4)) +
    labs(title = ttl, x = "Replicates, R", y = yl)

pb <- mk("b", "(b) detected in every replicate", expression(P(italic("possibly missing")~"|"~E)), 20)
pc <- mk("c", "(c) empty in every replicate", NULL, 20)

base <- theme_classic(base_size = 10) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.key.width = unit(18, "pt"),
        plot.title = element_text(size = 10, hjust = 0),
        axis.line = element_line(colour = "grey40", linewidth = 0.3),
        axis.ticks = element_line(colour = "grey40", linewidth = 0.3))

fig <- (pa + pb + pc) & base
fig <- fig + plot_layout(guides = "collect") & theme(legend.position = "bottom")

cat("kappa", kappa, "\n")
for (nu in c(NA, nus)) cat("nu", nu,
  " PM(n=R) R=1,3,5,10:", round(100*sapply(c(1,3,5,10), function(R) post(R,R,nu,PM)),1),
  " | PM(n=0) R=1,5,10,15,20:", round(100*sapply(c(1,5,10,15,20), function(R) post(R,0,nu,PM)),1), "\n")

ggsave("heterogeneity.pdf", fig, width = 9.5, height = 3.4)
ggsave("heterogeneity.png", fig, width = 9.5, height = 3.4, dpi = 300)
