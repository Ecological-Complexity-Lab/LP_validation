## ---------------------------------------------------------------------------
## SI figure: heterogeneous links. The realisation rate rho is drawn once per
## link from Beta(a, b) with mean rhobar and concentration nu = a + b, and held
## fixed across that link's replicates; the detectability 1 - eps_l stays common,
## so p_1 = rho(1 - eps_l) <= rho by construction.
##
## For a category with z_r = 1 the replicate factor is conditioned on the link
## being realised in at least one replicate, so both the numerator and the
## normaliser are expectations over the drawn rho (Eq. cumulative):
##
##   P(n | L_r = 1) = E_rho[ B*(n, R, rho u) ] / E_rho[ 1 - (1-rho)^R ],
##   B*(n,R,p) = C(R,n) p^n (1-p)^(R-n) - 1[n = 0] (1-rho)^R.
##
## Evidence: Y = 1, O_l = 0. Directional local axis (Eq. likelihood_dir).
## ---------------------------------------------------------------------------

library(ggplot2)
library(patchwork)

eY     <- 0.20
el     <- 0.30
f      <- 0.05
rhobar <- 0.15
u      <- 1 - el              # detectability, common to all links
nus    <- c(40, 10, 4)        # concentrations shown, plus the common rate
Rs     <- 1:20

## the eight categories as (z_Y, z_l, z_r)
sig <- rbind(c(1,1,1), c(1,1,0), c(1,0,1), c(1,0,0),
             c(0,1,1), c(0,1,0), c(0,0,1), c(0,0,0))
rownames(sig) <- c("recurrent","locally unique","possibly missing","phantom",
                   "model-elusive","weakly-supported","locally absent",
                   "possibly forbidden")
PM <- 3

## ---- replicate factor for z_r = 1, binomial coefficient dropped -------------
## Beta(a, b) with a < 1 is singular at 0, so integrate on a fine grid and check
## against adaptive quadrature (both agree to 1e-8 across the range shown).
gr <- seq(1e-9, 1 - 1e-9, length.out = 200001)

rep1 <- function(n, R, nu) {
  if (is.na(nu)) {                                   # common rate
    p <- rhobar * u
    num <- if (n == 0) (1 - p)^R - (1 - rhobar)^R else p^n * (1 - p)^(R - n)
    return(num / (1 - (1 - rhobar)^R))
  }
  w   <- dbeta(gr, rhobar * nu, (1 - rhobar) * nu)
  num <- if (n == 0) (1 - gr * u)^R - (1 - gr)^R else (gr * u)^n * (1 - gr * u)^(R - n)
  sum(w * num) / sum(w * (1 - (1 - gr)^R))
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
labs <- setNames(c("common rate", expression(nu == 40), expression(nu == 10),
                   expression(nu == 4)), lev)

## ---- (a) the distribution of the realisation rate across links -------------
gridp <- seq(1e-4, 0.6, length.out = 1200)
dens <- do.call(rbind, lapply(seq_along(nus), function(k) {
  data.frame(level = lev[k + 1], p = gridp,
             d = dbeta(gridp, rhobar * nus[k], (1 - rhobar) * nus[k]))
}))
dens$level <- factor(dens$level, levels = lev)

pa <- ggplot(dens, aes(p, d, colour = level)) +
  geom_vline(xintercept = rhobar, colour = cols[1], linewidth = 0.7) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = cols, labels = labs, drop = FALSE) +
  coord_cartesian(ylim = c(0, 11)) +
  labs(title = "(a) realisation rate across links",
       x = expression(rho), y = "density")

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
  scale_y_continuous(limits = c(0, 0.18), labels = pct) +
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
cat("\nE[rho | realised at least once]\n")
for (nu in nus) {
  w <- dbeta(gr, rhobar * nu, (1 - rhobar) * nu)
  for (R in c(1, 5, 20))
    cat(sprintf("  nu=%2d R=%2d  %.4f\n", nu, R,
        sum(w * gr * (1 - (1 - gr)^R)) / sum(w * (1 - (1 - gr)^R))))
}

ggsave("heterogeneity.pdf", fig, width = 7.2, height = 3.2)
ggsave("heterogeneity.png", fig, width = 7.2, height = 3.2, dpi = 300)
