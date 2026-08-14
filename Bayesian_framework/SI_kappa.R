## ---------------------------------------------------------------------------
## SI figure (fig:kappa): the maximum contextual confidence as a property of the
## model and the local method alone.
##   solid       symmetric rates, kappa = (1 - eps_Y)(1 - eps_l)      [sec:worked]
##   dashed      directional rates with O_l = 0 and false-detection rate f,
##               kappa = (1 - eps_Y)(1 - f) / ((1 - f) + eps_l)       [sec:grounding]
## The marked points are the two illustrative cases used in the text.
## ---------------------------------------------------------------------------

library(ggplot2)

eYs <- c(0.05, 0.20, 0.35)
f   <- 0.05
eLs <- seq(0, 0.6, by = 0.002)

kap_sym <- function(eY, eL) (1 - eY) * (1 - eL)
kap_dir <- function(eY, eL) (1 - eY) * (1 - f) / ((1 - f) + eL)

dat <- do.call(rbind, lapply(eYs, function(e) rbind(
  data.frame(eY = e, eL = eLs, kap = kap_sym(e, eLs), rates = "symmetric"),
  data.frame(eY = e, eL = eLs, kap = kap_dir(e, eLs), rates = "directional"))))
dat$eY    <- factor(dat$eY, levels = eYs)
dat$rates <- factor(dat$rates, levels = c("symmetric", "directional"))
dat$grp   <- interaction(dat$eY, dat$rates)

blues <- setNames(c("#0B3C6E", "#2A78D6", "#7FA9D4"), eYs)
labs  <- data.frame(eY = factor(eYs, levels = eYs), eL = 0.615,
                    kap = kap_dir(eYs, 0.6),
                    lab = paste0("epsilon[Y] == ", eYs))
pct <- function(x) paste0(round(100 * x), "%")

p <- ggplot(dat, aes(eL, kap, colour = eY, linetype = rates, group = grp)) +
  geom_line(linewidth = 0.8) +
  geom_text(data = labs, aes(eL, kap, label = lab), inherit.aes = FALSE,
            parse = TRUE, hjust = 0, vjust = 0.5, size = 3,
            colour = blues[as.character(eYs)]) +
  scale_colour_manual(values = blues, guide = "none") +
  scale_linetype_manual(values = c("solid", "22"), name = NULL,
                        labels = c(expression(epsilon[l] == f ~ "(symmetric)"),
                                   expression(f == 0.05 ~ "(directional)"))) +
  scale_y_continuous(limits = c(0, 1), labels = pct) +
  scale_x_continuous(limits = c(0, 0.78), breaks = seq(0, 0.6, 0.2)) +
  labs(x = expression("Local miss rate, " * epsilon[l]),
       y = expression("Maximum contextual confidence, " * kappa)) +
  theme_classic(base_size = 10) +
  theme(axis.line = element_line(colour = "grey40", linewidth = 0.3),
        axis.ticks = element_line(colour = "grey40", linewidth = 0.3),
        legend.position = "bottom")

cat(sprintf("eps_Y = 0.2, eps_l = 0.3: symmetric %.3f, directional %.3f\n",
            kap_sym(0.2, 0.3), kap_dir(0.2, 0.3)))
for (e in eYs)
  cat(sprintf("  eps_Y = %.2f: at eps_l = 0 both %.3f; at 0.6 symmetric %.3f, directional %.3f\n",
              e, kap_sym(e, 0), kap_sym(e, 0.6), kap_dir(e, 0.6)))

ggsave("SI_kappa.pdf", p, width = 5.4, height = 3.8)
ggsave("SI_kappa.png", p, width = 5.4, height = 3.8, dpi = 300)
