## ---------------------------------------------------------------------------
## SI figure: heterogeneous links. Under L_r = realisable in the replicates,
## rho enters the likelihood only through p_1 = rho(1 - eps_l), so it is enough
## to let p_1 vary. p_1 is drawn once per link from Beta(a, b) with mean pbar
## and concentration nu = a + b, and held fixed across that link's replicates,
## which makes the count of detections beta-binomial:
##
##   P(n | L_r = 1) = C(R,n) B(a+n, b+R-n) / B(a,b),   a = pbar*nu, b = (1-pbar)*nu
##
## Evidence: Y = 1, O_l = 0. Directional local axis (Eq. likelihood_dir).
## ---------------------------------------------------------------------------

library(ggplot2)
library(patchwork)

eY     <- 0.20
el     <- 0.30
f      <- 0.05
rho    <- 0.15
pbar   <- rho * (1 - el)      # 0.105, the common rate of Eq. cumulative
nus    <- c(100, 40, 10)      # concentrations shown, plus the common rate
Rs     <- 1:20

## the eight categories as (z_Y, z_l, z_r)
sig <- rbind(c(1,1,1), c(1,1,0), c(1,0,1), c(1,0,0),
             c(0,1,1), c(0,1,0), c(0,0,1), c(0,0,0))
rownames(sig) <- c("recurrent","locally unique","possibly missing","phantom",
                   "model-elusive","weakly-supported","locally absent",
                   "possibly forbidden")
PM <- 3

## ---- replicate factor for z_r = 1, binomial coefficient dropped -------------
rep1 <- function(n, R, nu) {
  if (is.na(nu)) pbar^n * (1 - pbar)^(R - n)                 # common rate
  else exp(lbeta(pbar * nu + n, (1 - pbar) * nu + R - n) -
           lbeta(pbar * nu, (1 - pbar) * nu))                # beta-binomial
}

post <- function(R, n, nu, which) {
  L <- numeric(8)
  for (i in 1:8) {
    lY <- if (sig[i,1] == 1) 1 - eY else eY
    lL <- if (sig[i,2] == 1) el else 1 - f           # O_l = 0, directional
    pr <- if (sig[i,3] == 1) rep1(n, R, nu) else f^n * (1 - f)^(R - n)
    L[i] <- lY * lL * pr
  }
  (L / sum(L))[which]
}

lev  <- c("common rate", paste0("nu == ", nus))
cols <- setNames(c("#4D4D4D", "#2A78D6", "#E07B39", "#7C5CBF"), lev)
labs <- setNames(c("common rate", expression(nu == 100), expression(nu == 40),
                   expression(nu == 10)), lev)

## ---- (a) the distribution of the detection rate across links ---------------
gridp <- seq(1e-4, 0.5, length.out = 1200)
dens <- do.call(rbind, lapply(seq_along(nus), function(k) {
  data.frame(level = lev[k + 1], p = gridp,
             d = dbeta(gridp, pbar * nus[k], (1 - pbar) * nus[k]))
}))
dens$level <- factor(dens$level, levels = lev)

pa <- ggplot(dens, aes(p, d, colour = level)) +
  geom_vline(xintercept = pbar, colour = cols[1], linewidth = 0.7) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = cols, labels = labs, drop = FALSE) +
  coord_cartesian(ylim = c(0, 16)) +
  labs(title = "(a) detection rate across links",
       x = expression(p[1]), y = "density")

## ---- (b) empty in every replicate ------------------------------------------
curves <- do.call(rbind, lapply(c(NA, nus), function(nu) {
  data.frame(level = if (is.na(nu)) lev[1] else paste0("nu == ", nu),
             R = Rs, p = sapply(Rs, function(R) post(R, 0, nu, PM)))
}))
curves$level <- factor(curves$level, levels = lev)

pct <- function(x) paste0(round(100 * x), "%")
pb <- ggplot(curves, aes(R, p, colour = level)) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = cols, labels = labs, drop = FALSE) +
  scale_y_continuous(limits = c(0, 0.35), labels = pct) +
  scale_x_continuous(breaks = scales::pretty_breaks(5)) +
  labs(title = "(b) empty in every replicate",
       x = "Replicates, R",
       y = expression(P(italic("possibly missing")~"|"~E)))

base <- theme_classic(base_size = 10) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.key.width = unit(18, "pt"),
        plot.title = element_text(size = 10, hjust = 0),
        axis.line = element_line(colour = "grey40", linewidth = 0.3),
        axis.ticks = element_line(colour = "grey40", linewidth = 0.3))

fig <- (pa + pb) & base
fig <- fig + plot_layout(guides = "collect") & theme(legend.position = "bottom")

## ---- reported values --------------------------------------------------------
cat("kappa =", (1 - eY) * (1 - f) / ((1 - f) + el), "\n")
for (nu in c(NA, nus)) {
  cat(sprintf("nu=%-7s", ifelse(is.na(nu), "common", as.character(nu))))
  for (R in c(1, 5, 10, 20)) cat(sprintf("  R%2d %.4f", R, post(R, 0, nu, PM)))
  cat("\n")
}
cat("\ndetections (n = R), possibly missing\n")
for (nu in c(NA, nus)) {
  cat(sprintf("nu=%-7s", ifelse(is.na(nu), "common", as.character(nu))))
  for (R in c(1, 5, 10, 20)) cat(sprintf("  R%2d %.4f", R, post(R, R, nu, PM)))
  cat("\n")
}
cat("\nphantom under empty replicates\n")
for (nu in c(NA, nus)) {
  cat(sprintf("nu=%-7s", ifelse(is.na(nu), "common", as.character(nu))))
  for (R in c(1, 5, 10, 20)) cat(sprintf("  R%2d %.4f", R, post(R, 0, nu, 4)))
  cat("\n")
}

ggsave("SI_heterogeneity.pdf", fig, width = 7.2, height = 3.2)
ggsave("SI_heterogeneity.png", fig, width = 7.2, height = 3.2, dpi = 300)
