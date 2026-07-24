# =============================================================================
# Box 2 - accumulation of confidence as within-system replicates accumulate.
# Two panels, one per consistent-evidence scenario:
#   (1) Detected in every replicate  -> k = R
#   (2) Empty in every replicate     -> k = 0
# All eight link categories, one distinct colour each. Observed within-network
# evidence (Yhat, O_local) is fixed by the `chosen` category.
#
# Model (Box 2, quantitative replicate axis):
#   score(c) = prior_c
#            * [1-eY if Yhat_c matches observed Yhat, else eY]
#            * [1-eL if Olocal_c matches observed Olocal, else eL]
#            * p_c^k * (1-p_c)^(R-k)
#   p_c = p1 = rho*(1-eL) if truly realised elsewhere (Orep=1); p_c = p0 = f if absent.
#   Posterior = score(c) / sum over all 8 categories.
# =============================================================================

library(ggplot2)

## ---- parameters -------------------------------------------------------------
eY   <- 0.20     # model error       (epsilon_Y)
eL   <- 0.30     # local-miss error  (epsilon_local)
rho  <- 0.15     # realisation rate
f    <- 0.05     # false-detection rate
Rmax <- 20       # max number of replicates on the x-axis
chosen <- "Possibly missing"   # sets the observed (Yhat, O_local)

## ---- the eight link categories ---------------------------------------------
## signature bits: Yhat (predicted), Olocal (observed locally), Orep (observed elsewhere)
cats <- data.frame(
  name   = c("Locally unique","Recurrent","Phantom","Possibly missing",
             "Possibly forbidden","Locally absent","Weak support","Model-elusive"),
  Yhat   = c(1,1,1,1, 0,0,0,0),
  Olocal = c(1,1,0,0, 0,0,1,1),
  Orep   = c(0,1,0,1, 0,1,0,1),
  stringsAsFactors = FALSE
)
cats$prior <- 1/nrow(cats)              # uniform prior
obs <- cats[cats$name == chosen, ]      # observed within-network evidence

## ---- posterior over the 8 categories for a given (R, k) ---------------------
posterior_at <- function(R, k) {
  p1 <- rho * (1 - eL)
  p0 <- f
  fY <- ifelse(obs$Yhat   == cats$Yhat,   1 - eY, eY)
  fL <- ifelse(obs$Olocal == cats$Olocal, 1 - eL, eL)
  p  <- ifelse(cats$Orep == 1, p1, p0)
  score <- cats$prior * fY * fL * p^k * (1 - p)^(R - k)
  score / sum(score)
}

## ---- build the data for both scenarios --------------------------------------
scenarios <- c(det = "Detected in every replicate (k = R)",
               emp = "Empty in every replicate (k = 0)")

df <- do.call(rbind, lapply(names(scenarios), function(sc) {
  do.call(rbind, lapply(0:Rmax, function(R) {
    k <- if (sc == "det") R else 0
    data.frame(scenario = scenarios[[sc]], R = R,
               category = cats$name, posterior = posterior_at(R, k),
               stringsAsFactors = FALSE)
  }))
}))

df$scenario <- factor(df$scenario, levels = scenarios)   # detected first
df$category <- factor(df$category, levels = cats$name)    # taxonomy order in legend

ceiling_val <- (1 - eY) * (1 - eL)   # common asymptote = (1-eY)(1-eL)

## ---- palette: one distinct colour per category (ColorBrewer Dark2) ----------
pal <- c(
  "Locally unique"     = "#1B9E77",
  "Recurrent"          = "#D95F02",
  "Phantom"            = "#7570B3",
  "Possibly missing"   = "#E7298A",
  "Possibly forbidden" = "#66A61E",
  "Locally absent"     = "#E6AB02",
  "Weak support"       = "#A6761D",
  "Model-elusive"      = "#666666"
)

## ---- plot -------------------------------------------------------------------
p <- ggplot(df, aes(R, posterior, colour = category, group = category)) +
  geom_hline(yintercept = ceiling_val, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 1) +
  geom_point(size = 1.4) +
  facet_wrap(~ scenario, ncol = 2) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(breaks = seq(0, Rmax, 5)) +
  labs(
    x = "number of replicates R (consistent evidence)",
    y = "posterior probability  P(C | E)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text       = element_text(face = "bold")
  )

print(p)

## ---- save -------------------------------------------------------------------
ggsave("box2_accumulation.pdf", p, width = 9.5, height = 4.2)
ggsave("box2_accumulation.png", p, width = 9.5, height = 4.2, dpi = 300)
