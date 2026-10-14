## ---------------------------------------------------------------------------
## SI figure (fig:SI_prominence): how a link's prominence elsewhere redirects
## accumulating replicate evidence. Replaces the earlier fig:SI_frequency,
## which swept the records at a fixed replicate count instead.
##
## THE LINK. The evidence is the one used throughout, Y = 1 and O_l = 0: the
## link is predicted and not recorded at the focal site. It is recorded in n of
## the R replicates, and n is the x-axis (Section frequency). Its ASSIGNED
## label is therefore possibly missing for every n >= 1, and phantom at n = 0.
##
## WHAT THE CURVES ARE. Not the assigned label, but the posterior that the
## link's TRUE label is one of the two the assignment could be hiding:
##   recurrent       the link is present at several sites, this one included,
##                   so the local non-detection was a miss
##   locally absent  the link is not realisable here and a model saying so
##                   would have been right, so there is nothing to find
## Both are read against the assigned label possibly missing, which holds the
## rest of the mass and is reported in the console.
##
## THE MECHANISM. Prominence enters in one place only: the link's own miss rate
## falls as its records accumulate (Eq. freq_miss, derived below). This is the
## only thing that couples the local axis to the replicate axis. With a common
## eps_l the posterior factorises over the three axes, so the replicate count
## alone says nothing about whether the non-detection was a miss, and both
## curves are flat in n; the console prints that comparison.
##   (a) weak where it occurs       one record in each recording replicate
##   (b) prominent where it occurs  ten records in each recording replicate
## The panels share their n = 0 end, where there are no records to read either
## way. The prior is uniform over the eight labels.
## ---------------------------------------------------------------------------

library(ggplot2)
library(patchwork)

## Figures are written to Bayesian_framework/bayesian_figures. Resolved here so
## the script works whether it is run from the project root or from inside
## Bayesian_framework/, and the folder is created if it is missing.
fig_dir <- if (dir.exists("Bayesian_framework"))
             file.path("Bayesian_framework", "bayesian_figures") else "bayesian_figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

eY <- 0.20; el <- 0.30; f <- 0.05; rho <- 0.15
R  <- 10                  # replicates in the system; n runs from 0 to R

p1 <- rho * (1 - el)      # a replicate records a realisable link, Eq. p1
b  <- el / (1 - el)       # fixes the average miss rate across links at eps_l

W_WEAK <- 1               # records per recording replicate, panel (a)
W_PROM <- 10              # records per recording replicate, panel (b)

sig <- rbind(c(1,1,1), c(1,1,0), c(1,0,1), c(1,0,0),
             c(0,1,1), c(0,1,0), c(0,0,1), c(0,0,0))
nm  <- c("recurrent", "locally unique", "possibly missing", "phantom",
         "model-elusive", "weakly supported", "locally absent",
         "possibly forbidden")

## ---- the link's own miss rate (Eq. freq_miss) ------------------------------
## THE MODEL. A link has its own interaction rate lambda: how often it happens
## where it is realisable. Give it an Exp(b) prior and let a replicate that can
## realise it hold Poisson(lambda) records. Then a sample misses the link when
## its Poisson draw is zero, with probability exp(-lambda), and averaging that
## over the prior gives b / (b + 1), which is eps_l exactly when b is set as
## above. So the common eps_l is this model with no records read; everything
## below is the same quantity after the records are read.
##
## THE ARGUMENTS. n replicates recorded the link, W records between them, out
## of R. Nothing else about the records matters: for Poisson counts n and W are
## sufficient, so ten records in one replicate and one in each of ten give the
## same answer only insofar as n and W agree.
##
## THE EMPTY REPLICATES. Each one is evidence too, and ambiguous: it either
## cannot realise the link (probability 1 - rho) or can and missed it
## (probability rho). j counts how many took the second reading, so the sum
## runs over every way of splitting the R - n empty replicates.
miss <- function(n, W) {
  j  <- 0:(R - n)
  ## P(j | the records), from the prior split times how well it explains them.
  ## The second factor is the Exp-Poisson evidence (b + n + j)^-(1 + W): a
  ## split that has the link sampled often (n + j large) yet yielding only W
  ## records is a split in which the link must be rare, so it is penalised.
  pj <- choose(R - n, j) * rho^j * (1 - rho)^(R - n - j) * (b + n + j)^-(1 + W)
  pj <- pj / sum(pj)
  ## Given the split, lambda is Gamma(1 + W, b + n + j), so the chance the
  ## focal site's sample also came up empty is E[exp(-lambda)], the bracket.
  ## It is near exp(-W / n), the records per realisable replicate: a link seen
  ## about once per sample keeps a miss rate near exp(-1) = 0.37 however many
  ## replicates record it, which is the ceiling panel (a) runs into, while ten
  ## records per replicate drive it to exp(-10) and panel (b) collapses.
  sum(pj * ((b + n + j) / (b + n + j + 1))^(1 + W))
}

## ---- posterior for Y = 1, O_l = 0 -------------------------------------------
## n of the R replicates recorded the link, w records in each of them. This is
## Eq. cumulative with one change: where the link is realisable in the
## replicates, its records set its miss rate.
## use_records = FALSE keeps the common eps_l, which is Eq. cumulative itself.
post <- function(n, w, use_records = TRUE) {
  W <- n * w
  L <- sapply(1:8, function(i) {
    fy <- if (sig[i,1] == 1) 1 - eY else eY
    ## the replicate records sharpen the miss rate only where they are real
    e  <- if (use_records && sig[i,3] == 1) miss(n, W) else el
    fl <- if (sig[i,2] == 1) e else 1 - f       # O_l = 0
    p  <- if (sig[i,3] == 1) p1 else f
    fy * fl * p^n * (1 - p)^(R - n)
  })
  L / sum(L)
}

## ---- the two curves at their ceiling ----------------------------------------
## Once n = R the labels that put the link nowhere are extinguished, the four
## realisable-elsewhere labels hold everything, and the miss rate e alone
## decides how they divide it. Writing them out shows what the figure is for:
## the same evidence splits between a miss and a true local absence according
## to e, and e is what the records set.
ceil_rec <- function(e) (1 - eY) * e / ((1 - f) + e)        # recurrent
ceil_abs <- function(e) eY * (1 - f) / ((1 - f) + e)        # locally absent

shown <- c(1, 7)                                 # recurrent, locally absent
## the Box 2 key, from CATEGORY_COLOUR in bayesian_empirical.R, so the curves
## carry the colours these categories have in the main map. Both lines are
## solid, so the key is all that separates them: the two colours share a hue
## and differ in lightness, which is what the map relies on as well.
pal <- setNames(c("#C05030", "#F9C8B0"), nm[shown])   # dark coral, pale peach

nvals <- 0:R
Y_TOP <- 0.30                                    # the curves never reach 1 - eps_Y

build <- function(w, ttl) {
  P <- sapply(nvals, function(n) post(n, w))
  d <- do.call(rbind, lapply(shown, function(i)
    data.frame(n = nvals, p = P[i, ], cat = nm[i])))
  d$cat <- factor(d$cat, levels = nm[shown])

  ## which label leads, named once per run of the x-axis it leads over. The
  ## changeover sits between two integers. Runs too short to hold a label are
  ## left unnamed rather than clipped, which is what happens in panel (b),
  ## where locally absent takes the lead at the first record.
  lead <- nm[shown][apply(P[shown, ], 2, which.max)]
  sw   <- which(diff(as.integer(factor(lead, levels = nm[shown]))) != 0)
  brk  <- (nvals[sw] + nvals[sw + 1]) / 2
  segs <- c(min(nvals), brk, max(nvals))
  mids <- (head(segs, -1) + tail(segs, -1)) / 2
  labs <- data.frame(x = mids, y = Y_TOP * 0.94, width = diff(segs),
                     lab = lead[sapply(mids, function(v) which.min(abs(nvals - v)))])
  labs <- labs[labs$width >= 1.5, ]

  ggplot(d, aes(n, p, colour = cat)) +
    geom_line(linewidth = 0.9) +
    geom_text(data = labs, aes(x, y, label = lab), inherit.aes = FALSE,
              size = 2.9, colour = "grey30", fontface = "italic") +
    scale_colour_manual(values = pal, name = NULL) +
    scale_y_continuous(limits = c(0, Y_TOP),
                       labels = function(x) paste0(round(100*x), "%")) +
    scale_x_continuous(breaks = seq(0, R, 2)) +
    labs(title = ttl,
         x = expression("Replicates recording the link, " * italic(n)),
         y = expression(P("true label" ~ "|" ~ E))) +
    theme_classic(base_size = 10) +
    ## two-line title at 9 pt: the record level has to be stated, and on one
    ## line at 10 pt it runs past the panel and is clipped
    theme(plot.title = element_text(size = 10, hjust = 0, lineheight = 1.15),
          axis.line = element_line(colour = "grey40", linewidth = 0.3),
          axis.ticks = element_line(colour = "grey40", linewidth = 0.3),
          legend.position = "bottom")
}

## plain strings, not bquote: there is no maths in the titles, and plotmath
## neither wraps nor measures as compactly
pa <- build(W_WEAK, sprintf("(a) weak where it occurs",
                            W_WEAK))
pb <- build(W_PROM, sprintf("(b) prominent where it occurs",
                            W_PROM)) +
  labs(y = NULL)

fig <- pa + pb + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

## ---- reported values --------------------------------------------------------
cat(sprintf("\nR = %d, eps_Y = %.2f, eps_l = %.2f, f = %.2f, rho = %.2f\n",
            R, eY, el, f, rho))
report <- function(w, lbl) {
  cat(sprintf("\n%s\n", lbl))
  cat(sprintf("%4s %9s %11s %14s %17s %9s  leads\n", "n", "miss rate",
              nm[1], nm[7], nm[3], nm[4]))
  for (n in nvals) {
    p <- post(n, w)
    cat(sprintf("%4d %9.4f %11.3f %14.3f %17.3f %9.3f  %s\n", n, miss(n, n * w),
                p[1], p[7], p[3], p[4], nm[shown][which.max(p[shown])]))
  }
  P <- sapply(nvals, function(n) post(n, w))
  e <- miss(R, R * w)
  cat(sprintf("  %s peaks at %.3f, n = %d\n", nm[1],
              max(P[1, ]), nvals[which.max(P[1, ])]))
  cat(sprintf("  at n = R: miss rate %.4f against exp(-w) = %.4f, %s %.3f, %s %.3f\n",
              e, exp(-w), nm[1], ceil_rec(e), nm[7], ceil_abs(e)))
}
report(W_WEAK, "(a) one record per recording replicate")
report(W_PROM, "(b) ten records per recording replicate")

## the same posteriors with the common eps_l, flat in n because the three axes
## then factorise: this is what the records buy us
cat("\nEq. cumulative, common eps_l\n")
cat(sprintf("%4s %11s %14s %17s\n", "n", nm[1], nm[7], nm[3]))
for (n in nvals) {
  p <- post(n, W_WEAK, use_records = FALSE)
  cat(sprintf("%4d %11.3f %14.3f %17.3f\n", n, p[1], p[7], p[3]))
}
cat(sprintf("  implied by eps_l = %.2f: %s %.3f, %s %.3f\n",
            el, nm[1], ceil_rec(el), nm[7], ceil_abs(el)))

## the two curves are a fixed ratio of each other, since they differ only in
## whether the model is right, so the figure also fixes possibly missing
cat(sprintf("\n%s / %s and the label the link is assigned\n", nm[7], nm[1]))
for (w in c(W_WEAK, W_PROM)) {
  e <- miss(R, R * w)
  cat(sprintf("  w = %2d: e = %.4f, ratio %.3f, %s %.3f\n", w, e,
              ceil_abs(e) / ceil_rec(e), nm[3], (1 - eY) * (1 - f) / ((1 - f) + e)))
}

## ---- check Eq. freq_miss against direct numerical integration --------------
## Carried over from the retired fig:SI_frequency script: the closed form above
## is the only non-obvious step, so it is verified against the integral it came
## from, E[exp(-lambda) | records], over the eight record patterns the figure
## relies on. w is the full vector of records, one entry per replicate.
##
## Written in logs, rescaled by its peak and integrated over a finite window
## rather than (0, Inf). Both changes are needed at the record levels panel (b)
## uses: with ten records in each of ten replicates the integrand's total mass
## is around 1e-20, and integrate() compares that against its default absolute
## tolerance of about 1e-10, declares convergence after two subdivisions and
## returns a value with an error bar larger than itself.
miss_num <- function(w) {
  n <- sum(w > 0); W <- sum(w)
  ## log integrand: Exp(b) prior, Poisson records where the link is realisable,
  ## and the two readings of an empty replicate. k = 1 inserts the exp(-lambda)
  ## whose posterior mean is the miss rate; k = 0 is the normalising constant.
  lg <- function(x, k) log(b) - (b + k) * x +
    sum(ifelse(w == 0, log((1 - rho) + rho * exp(-x)),
               log(rho) + dpois(w, x, log = TRUE)))
  peak <- max(W / max(n, 1), 0.5)            # the rate the records point to
  hi   <- peak + 40 * sqrt(peak + 1) + 40    # past it by far more than needed
  M    <- lg(peak, 0)                        # rescale so the peak is 1
  g    <- function(x, k) sapply(x, function(v) exp(lg(v, k) - M))
  integrate(g, 0, hi, k = 1, rel.tol = 1e-12, subdivisions = 2000)$value /
    integrate(g, 0, hi, k = 0, rel.tol = 1e-12, subdivisions = 2000)$value
}
for (w in list(rep(0, R), c(1, rep(0, R - 1)), c(10, rep(0, R - 1)),
               c(5, 5, rep(0, R - 2)), c(2, 1, 1, rep(0, R - 3)),
               rep(W_WEAK, R), rep(W_PROM, R), c(rep(W_PROM, R - 1), 0))) {
  ## one entry per replicate, empties included: a short vector would describe a
  ## smaller system than the miss() being checked against it
  stopifnot(length(w) == R)
  stopifnot(abs(miss(sum(w > 0), sum(w)) - miss_num(w)) < 1e-9)
}

fig
ggsave(file.path(fig_dir, "SI_prominence.pdf"), fig, width = 7.6, height = 3.6)
ggsave(file.path(fig_dir, "SI_prominence.png"), fig, width = 7.6, height = 3.6, dpi = 300)
