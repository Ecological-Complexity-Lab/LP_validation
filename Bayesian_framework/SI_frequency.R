## ---------------------------------------------------------------------------
## SI figure (fig:SI_frequency): what interaction frequencies add to a fixed
## piece of evidence. The evidence is the one used throughout: the link is
## predicted (Y = 1) and not recorded locally (w_l = 0). It is recorded in all
## five replicates, and what varies is how many records each replicate holds
## (Section frequency).
## Frequencies enter in two places:
##   1. a false record of w records needs w separate errors, f^w
##   2. the link's own miss rate falls as its records accumulate
## Two categories are drawn, the two in which the link is unrealised at the
## focal site and realisable in the replicates: possibly missing (the model is
## right) and locally absent (the model is wrong). Together they approach 100%,
## split by the model error alone.
## ---------------------------------------------------------------------------

library(ggplot2)

## Figures are written to Bayesian_framework/bayesian_figures. Resolved here so
## the script works whether it is run from the project root or from inside
## Bayesian_framework/, and the folder is created if it is missing.
fig_dir <- if (dir.exists("Bayesian_framework"))
             file.path("Bayesian_framework", "bayesian_figures") else "bayesian_figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

eY <- 0.20; el <- 0.30; f <- 0.05; rho <- 0.15
R  <- 5

p1 <- rho * (1 - el)      # a replicate records a realisable link, Eq. p1
b  <- el / (1 - el)       # fixes the average miss rate across links at eps_l

sig <- rbind(c(1,1,1), c(1,1,0), c(1,0,1), c(1,0,0),
             c(0,1,1), c(0,1,0), c(0,0,1), c(0,0,0))
nm  <- c("recurrent", "locally unique", "possibly missing", "phantom",
         "model-elusive", "weakly supported", "locally absent",
         "possibly forbidden")

## ---- the link's own miss rate ----------------------------------------------
## n replicates recorded the link, with W records between them. With no
## records this returns the common eps_l.
miss <- function(n, W) ((b + n) / (b + n + 1))^(1 + W)

## ---- posterior with record counts ------------------------------------------
## wl = records at the focal site, w = records in each replicate, Y = 1.
## q is optional: the probability that a realised link recorded at a site is
## recorded w times. Left out, a recording replicate counts as recorded,
## whatever its count (Eq. frequency). Supplied, the count in each replicate
## enters as well (Eq. frequency_full).
post_freq <- function(wl, w, q = NULL) {
  n <- sum(w > 0); W <- sum(w)
  qw <- if (is.null(q)) 1 else prod(q(w[w > 0]))
  L <- sapply(1:8, function(i) {
    fy <- if (sig[i,1] == 1) 1 - eY else eY
    ## the replicate records sharpen the miss rate only where they are real
    e  <- if (sig[i,3] == 1) miss(n, W) else el
    fl <- if (sig[i,2] == 1) { if (wl == 0) e else (1 - e) * (if (is.null(q)) 1 else q(wl)) }
          else               { if (wl == 0) 1 - f else f^wl }
    fr <- if (sig[i,3] == 1) p1^n * qw * (1 - p1)^(R - n)
          else               f^W * (1 - f)^(R - n)
    fy * fl * fr
  })
  L / sum(L)
}

shown <- c(3, 7)                                 # unrealised here, realisable elsewhere
## locally absent is not in the palette of the other SI figures; replace its
## colour with the one this category has in the main-text figures
pal <- setNames(c("#E07B39", "slateblue3"), nm[shown])

kappa <- (1 - eY) * (1 - f) / ((1 - f) + el)     # ceiling without frequencies

x <- 0:10                                        # records in each replicate
all_rep <- function(v) rep(v, R)                 # the same records in every replicate

P <- sapply(x, function(v) post_freq(0, all_rep(v)))
d <- do.call(rbind, lapply(shown, function(i)
  data.frame(w = x, p = P[i, ], cat = nm[i])))
d$cat <- factor(d$cat, levels = nm[shown])

refs <- data.frame(y = c(1 - eY, kappa, eY),
                   lab = c("1 - epsilon[Y]", "kappa", "epsilon[Y]"))

fig <- ggplot(d, aes(w, p, colour = cat)) +
  geom_hline(data = refs, aes(yintercept = y), linetype = "dotted",
             colour = "grey60", linewidth = 0.35) +
  geom_text(data = refs, aes(x = max(x), y = y + 0.035, label = lab),
            inherit.aes = FALSE, parse = TRUE, hjust = 1,
            size = 2.9, colour = "grey30") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.1) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(round(100*x), "%")) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(
       x = "Records in each replicate", y = expression(P(C ~ "|" ~ E))) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0),
        axis.line = element_line(colour = "grey40", linewidth = 0.3),
        axis.ticks = element_line(colour = "grey40", linewidth = 0.3),
        legend.position = "bottom")

## ---- reported values --------------------------------------------------------
cat(sprintf("\nkappa = %.3f, ceiling with frequencies = %.3f\n", kappa, 1 - eY))
report <- function(lbl, w_of) {
  cat(sprintf("\n%s\n", lbl))
  cat(sprintf("%8s %9s %17s %15s %10s %8s\n", "records", "miss rate",
              nm[3], nm[7], nm[1], nm[4]))
  for (v in c(0, 1, 2, 5, 10)) {
    w <- w_of(v); p <- post_freq(0, w)
    cat(sprintf("%8d %9.3f %17.3f %15.3f %10.3f %8.3f\n", v,
                miss(sum(w > 0), sum(w)), p[3], p[7], p[1], p[4]))
  }
}
report("all five replicates record the link (the figure)", all_rep)
report("one replicate of five records the link (text only)",
       function(v) c(v, rep(0, R - 1)))

## records at the focal site, none in the replicates: how likely the local
## record is false (the four categories with z_l = 0)
cat("\nlocal records, no replicate records: P(local record is false)\n")
no_local <- which(sig[, 2] == 0)
for (wl in 1:2)
  cat(sprintf("  %d local record%s: %.4f\n", wl, if (wl > 1) "s" else " ",
              sum(post_freq(wl, rep(0, R))[no_local])))

## the optional count term, under the count distribution that the exponential
## rate implies, q(w) = eps_l (1 - eps_l)^(w - 1): single records then count
## against the link
q_exp <- function(w) el * (1 - el)^(w - 1)
cat("\none record in each of five replicates: phantom\n")
cat(sprintf("  without the count term: %.3f\n", post_freq(0, all_rep(1))[4]))
cat(sprintf("  with q from the exponential: %.3f\n", post_freq(0, all_rep(1), q_exp)[4]))
cat("ten records in two replicates, with q: possibly missing\n")
cat(sprintf("  5 and 5: %.3f   9 and 1: %.3f\n",
            post_freq(0, c(5, 5, 0, 0, 0), q_exp)[3], post_freq(0, c(9, 1, 0, 0, 0), q_exp)[3]))


print(fig)
ggsave(file.path(fig_dir, "SI_frequency.pdf"), fig, width = 4.2, height = 3.5)
ggsave(file.path(fig_dir, "SI_frequency.png"), fig, width = 4.2, height = 3.5, dpi = 300)
