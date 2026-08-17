## ---------------------------------------------------------------------------
## Main-text figure (fig:main_accumulation): what replicate evidence buys,
## and what it cannot.
## (a) Confidence accumulates with consistent detections but stops at kappa,
##     while the aggregate feasibility confidence phi is not bounded by it.
## (b) kappa = (1 - eps_Y)(1 - eps_l) is a property of the model and the local
##     method alone; the rule in (a) is one point on these curves.
## Symmetric rates throughout, matching the likelihood printed in the section.
## ---------------------------------------------------------------------------

library(ggplot2)
library(patchwork)

## Figures are written to Bayesian_framework/bayesian_figures. Resolved here so
## the script works whether it is run from the project root or from inside
## Bayesian_framework/, and the folder is created if it is missing.
fig_dir <- if (dir.exists("Bayesian_framework"))
             file.path("Bayesian_framework", "bayesian_figures") else "bayesian_figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

eY <- 0.20; eL <- 0.30            # model and local error rates
p1 <- 0.105; p0 <- 0.05           # per-replicate recording probability
Rs <- 0:20

## eight categories as (z_Y, z_l, z_r); evidence is Y = 1, O_l = 0
sig <- rbind(c(1,1,1), c(1,1,0), c(1,0,1), c(1,0,0),
             c(0,1,1), c(0,1,0), c(0,0,1), c(0,0,0))
rownames(sig) <- c("recurrent","locally unique","possibly missing","phantom",
                   "model-elusive","weakly-supported","locally absent",
                   "possibly forbidden")

## model and local factors for the evidence Y = 1, O_l = 0 (symmetric local axis)
w <- ifelse(sig[,1] == 1, 1 - eY, eY) * ifelse(sig[,2] == 0, 1 - eL, eL)

post <- function(R) {                       # detected in every replicate
  if (R == 0) return(w / sum(w))            # no replicate term at all
  ## z_r = 1 means the link is realisable in the replicates, so the count is
  ## plainly binomial and no conditioning is applied (Eqs. erep, cumulative)
  pr <- ifelse(sig[,3] == 1, p1^R, p0^R)
  L  <- w * pr
  L / sum(L)
}

P   <- t(sapply(Rs, post)); colnames(P) <- rownames(sig)
kap <- (1 - eY) * (1 - eL)                  # 0.56

## two pairs of categories, each pair differing only in the replicate bit
lev <- c("Possibly missing","Phantom","Recurrent","Locally unique")
key <- c("Possibly missing (1, 0, 1)", "Phantom (1, 0, 0)",
         "Recurrent (1, 1, 1)",        "Locally unique (1, 1, 0)")
dat_a <- rbind(
  data.frame(R = Rs, y = P[,"possibly missing"],   k = lev[1]),
  data.frame(R = Rs, y = P[,"phantom"],            k = lev[2]),
  data.frame(R = Rs, y = P[,"recurrent"],          k = lev[3]),
  data.frame(R = Rs, y = P[,"locally unique"],     k = lev[4]))
dat_a$k <- factor(dat_a$k, levels = lev)
cols <- setNames(c("#E07B39","#E07B39","#7C5CBF","#7C5CBF"), lev)
ltys <- setNames(c("solid","22","solid","22"), lev)

phi <- data.frame(R = Rs, y = 1 - P[,"phantom"] - P[,"possibly forbidden"])

pct <- function(x) paste0(round(100 * x), "%")

pa <- ggplot(dat_a, aes(R, y, colour = k, linetype = k)) +
  geom_hline(yintercept = kap, linetype = "22", colour = "grey45",
             linewidth = 0.4) +
  annotate("text", x = 0, y = kap, label = "kappa == 0.56", parse = TRUE,
           hjust = 0, vjust = -0.6, size = 3, colour = "grey35",
           fontface = "bold") +
  geom_line(data = phi, aes(R, y), inherit.aes = FALSE, colour = "#2A9D8F",
            linewidth = 0.8, linetype = "12") +
  annotate("text", x = 20, y = 1, label = "feasibility confidence", hjust = 1,
           vjust = 1.6, size = 2.9, colour = "#2A9D8F", fontface = "bold") +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols, labels = key,
                      name = expression(bold("Category ") *
                        bold("(") * bolditalic(Y) * bold(", ") *
                        bolditalic(O)[bold(l)] * bold(", ") *
                        bolditalic(O)[bold(r)] * bold(")"))) +
  scale_linetype_manual(values = ltys, labels = key,
                        name = expression(bold("Category ") *
                          bold("(") * bolditalic(Y) * bold(", ") *
                          bolditalic(O)[bold(l)] * bold(", ") *
                          bolditalic(O)[bold(r)] * bold(")"))) +
  scale_y_continuous(limits = c(0, 1), labels = pct) +
  scale_x_continuous(breaks = seq(0, 20, 5)) +
  labs(title = "evidence accumulates to a maximum",
       x = "Replicates recording the link, R", y = "Posterior probability") +
  guides(colour = guide_legend(nrow = 2, title.position = "top"),
         linetype = guide_legend(nrow = 2, title.position = "top"))

## ---- (b) the ceiling itself -------------------------------------------------
## ---- (b) what an informative prior does to the same evidence ---------------
## Binary replicate evidence (O_r = 1 after R = 5 replicates), as in the
## framework as first presented; both priors rise with the pair's degree, the
## local one held just below the regional one.
Rb   <- 5
eP_b <- (1 - p1)^Rb              # eps_r^+(R), a realisable link recorded in none
eM_b <- 1 - (1 - p0)^Rb          # eps_r^-(R), an unrealisable one recorded in some
kfrac <- 0.9                     # pi_l = kfrac * pi_r

post_pi <- function(pr) {
  pl  <- kfrac * pr
  pri <- sapply(1:8, function(i)
    prod(c(0.5, pl, pr)^sig[i, ] * (1 - c(0.5, pl, pr))^(1 - sig[i, ])))
  L <- w * ifelse(sig[, 3] == 1, 1 - eP_b, eM_b) * pri
  L / sum(L)
}

prs <- seq(0.02, 0.95, length.out = 400)
Pb  <- sapply(prs, post_pi)
rownames(Pb) <- rownames(sig)
dat_b <- rbind(
  data.frame(pr = prs, y = Pb["possibly missing", ], k = lev[1]),
  data.frame(pr = prs, y = Pb["phantom", ],         k = lev[2]),
  data.frame(pr = prs, y = Pb["recurrent", ],       k = lev[3]),
  data.frame(pr = prs, y = Pb["locally unique", ],  k = lev[4]))
dat_b$k <- factor(dat_b$k, levels = lev)

## where the leading category changes
ld  <- rownames(sig)[apply(Pb[1:4, ], 2, which.max)]
brk <- prs[which(diff(as.integer(factor(ld))) != 0)]
seg <- c(min(prs), brk, max(prs))
labs_b <- data.frame(x = (head(seg, -1) + tail(seg, -1)) / 2, y = 0.95,
                     lab = ld[sapply((head(seg, -1) + tail(seg, -1)) / 2,
                                     function(v) which.min(abs(prs - v)))])

pb <- ggplot(dat_b, aes(pr, y, colour = k, linetype = k)) +
  geom_vline(xintercept = brk, linetype = "dotted", colour = "grey60",
             linewidth = 0.35) +
  geom_line(linewidth = 0.8) +
  geom_text(data = labs_b, aes(x, y, label = lab), inherit.aes = FALSE,
            size = 2.9, colour = "grey30", fontface = "italic") +
  scale_colour_manual(values = cols, guide = "none") +
  scale_linetype_manual(values = ltys, guide = "none") +
  scale_y_continuous(limits = c(0, 1), labels = pct) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(title = "an informative prior changes which error explains it",
       x = expression("Prior that the link is realisable in the replicates, " * pi[r]),
       y = "Posterior probability")

base <- theme_classic(base_size = 10) +
  theme(legend.position = "bottom", legend.key.width = unit(20, "pt"),
        legend.text = element_text(face = "bold", size = 7.5),
        legend.title = element_text(size = 8.5),
        plot.title = element_text(size = 9, hjust = 0, face = "bold"),
        ## panel letter: a bare "a" / "b", set larger than the title. Carried by
        ## patchwork's tag rather than the title string, so the letter and the
        ## title can be sized independently.
        plot.tag = element_text(size = 13, face = "bold", hjust = 0, vjust = 1),
        plot.tag.position = c(0, 1),
        axis.line = element_line(colour = "grey40", linewidth = 0.3),
        axis.ticks = element_line(colour = "grey40", linewidth = 0.3))

fig <- ((pa + pb) & base) + plot_annotation(tag_levels = "a")
ggsave(file.path(fig_dir, "main_accumulation.pdf"), fig, width = 9.0, height = 3.8)
ggsave(file.path(fig_dir, "main_accumulation.png"), fig, width = 9.0, height = 3.8, dpi = 300)

cat(sprintf("kappa = %.3f\n", kap))
print(round(100 * data.frame(R = Rs, pm = P[,"possibly missing"],
      ph = P[,"phantom"], phi = 1 - P[,"phantom"] - P[,"possibly forbidden"])[c(1,2,4,6,11,21),-1], 1))
