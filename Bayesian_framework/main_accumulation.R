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

eY <- 0.20; eL <- 0.30            # model and local error rates
p1 <- 0.105; p0 <- 0.05           # per-replicate recording probability
Rs <- 0:20

## eight categories as (z_Y, z_l, z_r); evidence is Y = 1, O_l = 0
sig <- rbind(c(1,1,1), c(1,1,0), c(1,0,1), c(1,0,0),
             c(0,1,1), c(0,1,0), c(0,0,1), c(0,0,0))
rownames(sig) <- c("recurrent","locally unique","possibly missing","phantom",
                   "model-elusive","weakly-supported","locally absent",
                   "possibly forbidden")

post <- function(R) {                       # detected in every replicate
  w  <- ifelse(sig[,1] == 1, 1 - eY, eY) *
        ifelse(sig[,2] == 0, 1 - eL, eL)    # symmetric local axis
  pr <- ifelse(sig[,3] == 1, p1, p0)^R
  L  <- w * pr
  L / sum(L)
}

P   <- t(sapply(Rs, post)); colnames(P) <- rownames(sig)
kap <- (1 - eY) * (1 - eL)                  # 0.56

## two pairs of categories, each pair differing only in the replicate bit
lev <- c("Possibly missing","Phantom","Locally absent","Possibly forbidden")
key <- c("Possibly missing (1, 0, 1)", "Phantom (1, 0, 0)",
         "Locally absent (0, 0, 1)",  "Possibly forbidden (0, 0, 0)")
dat_a <- rbind(
  data.frame(R = Rs, y = P[,"possibly missing"],   k = lev[1]),
  data.frame(R = Rs, y = P[,"phantom"],            k = lev[2]),
  data.frame(R = Rs, y = P[,"locally absent"],     k = lev[3]),
  data.frame(R = Rs, y = P[,"possibly forbidden"], k = lev[4]))
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
  labs(title = "(a) evidence accumulates to a maximum",
       x = "Replicates recording the link, R", y = "Posterior probability") +
  guides(colour = guide_legend(nrow = 2, title.position = "top"),
         linetype = guide_legend(nrow = 2, title.position = "top"))

## ---- (b) the ceiling itself -------------------------------------------------
eYs <- c(0.05, 0.20, 0.35)
dat_b <- do.call(rbind, lapply(eYs, function(e)
  data.frame(eY = e, eL = seq(0, 0.6, by = 0.005),
             kap = (1 - e) * (1 - seq(0, 0.6, by = 0.005)))))
dat_b$eY <- factor(dat_b$eY, levels = eYs)
labs_b <- data.frame(eY = factor(eYs, levels = eYs), eL = 0.015,
                     kap = (1 - eYs),
                     lab = paste0("epsilon[Y] == ", eYs))
blues <- setNames(c("#0B3C6E","#2A78D6","#7FA9D4"), eYs)

pb <- ggplot(dat_b, aes(eL, kap, colour = eY)) +
  geom_line(linewidth = 0.8) +
  geom_text(data = labs_b, aes(label = lab), parse = TRUE, hjust = 0,
            vjust = -0.8, size = 2.9, show.legend = FALSE, fontface = "bold") +
  annotate("point", x = eL, y = kap, size = 2, colour = "grey25") +
  annotate("text", x = eL, y = kap, label = "kappa~plain(\"in (a)\")", parse = TRUE, fontface = "bold", hjust = 0.5,
           vjust = 1.9, size = 2.9, colour = "grey25") +
  scale_colour_manual(values = blues, guide = "none") +
  scale_y_continuous(limits = c(0, 1.03), labels = pct) +
  labs(title = "(b) the maximum is set by the model and the local method",
       x = expression("Local miss rate, " * epsilon[l]),
       y = expression("Maximum contextual confidence, " * kappa))

base <- theme_classic(base_size = 10) +
  theme(legend.position = "bottom", legend.key.width = unit(20, "pt"),
        legend.text = element_text(face = "bold", size = 7.5),
        legend.title = element_text(size = 8.5),
        plot.title = element_text(size = 9, hjust = 0, face = "bold"),
        axis.line = element_line(colour = "grey40", linewidth = 0.3),
        axis.ticks = element_line(colour = "grey40", linewidth = 0.3))

fig <- (pa + pb) & base
ggsave("main_accumulation.pdf", fig, width = 9.0, height = 3.8)
ggsave("main_accumulation.png", fig, width = 9.0, height = 3.8, dpi = 300)

cat(sprintf("kappa = %.3f\n", kap))
print(round(100 * data.frame(R = Rs, pm = P[,"possibly missing"],
      ph = P[,"phantom"], phi = 1 - P[,"phantom"] - P[,"possibly forbidden"])[c(1,2,4,6,11,21),-1], 1))
