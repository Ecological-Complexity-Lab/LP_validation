# =============================================================================
# Box 2 - accumulation of confidence as within-system replicates accumulate.
# Two panels, one per consistent-evidence scenario:
#   (a) detected in every replicate  -> n = R
#   (b) empty in every replicate     -> n = 0
#
# Within each panel two summaries of the same replicate evidence are compared:
#   count : the full count n out of R   -> p_c^n (1-p_c)^(R-n)
#   bit   : O_rep = 1[n >= 1]           -> 1-(1-p_c)^R if O_rep = 1,
#                                          (1-p_c)^R   otherwise
# Under n = 0 the two coincide exactly, since (1-p)^R is both, so the dashed
# curve lies on top of the solid one in panel (b).
#
# Error rates are directional, as in "Grounding the error rates in sampling
# and realisation": the local axis errs at eL on a realised link and at f on
# an absent one; the model axis is symmetric at eY.
#   p_c = p1 = rho*(1-eL) if realised elsewhere (Orep = 1); p_c = p0 = f if not.
# =============================================================================

library(ggplot2)

## ---- parameters -------------------------------------------------------------
eY   <- 0.20     # model error            (epsilon_Y, symmetric)
eL   <- 0.30     # local miss rate        (epsilon_l)
f    <- 0.05     # false-detection rate   (f)
rho  <- 0.15     # realisation rate       (rho)
Rmax <- 30       # max number of replicates on the x-axis
chosen <- "Possibly missing"   # sets the observed (Yhat, O_local)

## ---- the eight link categories ---------------------------------------------
## signature bits: Yhat (predicted), Olocal (observed locally), Orep (elsewhere)
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

p1 <- rho * (1 - eL)
p0 <- f

## ---- within-network factors (model symmetric, local directional) -----------
fY <- ifelse(obs$Yhat == cats$Yhat, 1 - eY, eY)
fL <- ifelse(cats$Olocal == 1,
             ifelse(obs$Olocal == 1, 1 - eL, eL),    # link realised locally
             ifelse(obs$Olocal == 0, 1 - f,  f))     # link absent locally
w  <- cats$prior * fY * fL
p  <- ifelse(cats$Orep == 1, p1, p0)

## ---- replicate factor -------------------------------------------------------
## For categories with Orep = 1 the rates are conditional on the link being
## realised in at least one replicate: divide by Z = 1-(1-rho)^R and, at n = 0,
## remove the runs in which it was realised nowhere (Eqs. erep and cumulative).
Zc  <- function(R) ifelse(cats$Orep == 1, 1 - (1 - rho)^R, 1)
corr <- function(n, R) ifelse(cats$Orep == 1 & n == 0, (1 - rho)^R, 0)

eps_present <- function(R) ((1 - p1)^R - (1 - rho)^R) / (1 - (1 - rho)^R)
eps_absent  <- function(R) 1 - (1 - p0)^R

posterior_at <- function(R, n, summary = c("count","bit")) {
  summary <- match.arg(summary)
  ## R = 0: no replicate term at all, so every category keeps only w. The eight
  ## categories pair up by Orep, so each pair splits evenly and the shown pair
  ## starts at kappa/2.
  if (R == 0) return(w / sum(w))
  rep_lik <- if (summary == "count") {
    (p^n * (1 - p)^(R - n) - corr(n, R)) / Zc(R)   # binom coeff cancels
  } else {
    eP <- eps_present(R); eM <- eps_absent(R)
    if (n >= 1) ifelse(cats$Orep == 1, 1 - eP, eM)
      else      ifelse(cats$Orep == 1, eP, 1 - eM)
  }
  score <- w * rep_lik
  score / sum(score)
}

## ---- build the data ---------------------------------------------------------
scenarios <- c(det = "(a) detected in every replicate (n = R)",
               emp = "(b) empty in every replicate (n = 0)")
summaries <- c(count = "count", bit = "bit")
sum_labs  <- c(expression("full count " * italic(n)),
               expression("single bit " * bold(1)*"["*italic(n) >= 1*"]"))

df <- do.call(rbind, lapply(names(scenarios), function(sc) {
  do.call(rbind, lapply(names(summaries), function(sm) {
    do.call(rbind, lapply(0:Rmax, function(R) {
      n <- if (sc == "det") R else 0
      data.frame(scenario = scenarios[[sc]], summary = summaries[[sm]], R = R,
                 category = cats$name, posterior = posterior_at(R, n, sm),
                 stringsAsFactors = FALSE)
    }))
  }))
}))

shown <- c("Possibly missing","Phantom")
df <- df[df$category %in% shown, ]
df$category <- factor(df$category, levels = shown)
df$summary  <- factor(df$summary,  levels = summaries)
df$scenario <- factor(df$scenario, levels = scenarios)

kappa <- (1 - eY) * (1 - f) / ((1 - f) + eL)   # directional ceiling

## annotation for panel (b): the two summaries are the same number at n = 0
ann <- data.frame(scenario = factor(scenarios[["emp"]], levels = scenarios),
                  R = Rmax / 2, posterior = 0.92,
                  label = "the two summaries coincide",
                  stringsAsFactors = FALSE)

## ---- reported values --------------------------------------------------------
cat(sprintf("kappa = %.4f\n", kappa))
key <- df[df$R %in% c(1,3,5,10,20,26), ]
print(key[order(key$scenario, key$summary, key$category, key$R), ], row.names = FALSE)

## ---- palette (matched to the published figure) ------------------------------
pal <- c("Possibly missing" = "#E07B39", "Phantom" = "#7C5CBF")

## ---- plot -------------------------------------------------------------------
p_fig <- ggplot(df, aes(R, posterior, colour = category, linetype = summary,
                        shape = summary, group = interaction(category, summary))) +
  geom_hline(yintercept = kappa, linetype = "dashed", colour = "grey55",
             linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.3, stroke = 0.5, fill = "white") +
  geom_text(data = ann, aes(R, posterior, label = label), inherit.aes = FALSE,
            size = 3.1, colour = "grey30", hjust = 0.5) +
  facet_wrap(~ scenario, ncol = 2) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_linetype_manual(values = c("solid","22"), name = "replicate evidence",
                        labels = sum_labs) +
  scale_shape_manual(values = c(16, 21), name = "replicate evidence",
                     labels = sum_labs) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(breaks = seq(0, Rmax, 10)) +
  labs(x = "number of replicates R (consistent evidence)",
       y = "posterior probability  P(C | E)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "right",
        legend.key.width = unit(1.1, "cm"),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text       = element_text(face = "bold", size = 10))

## ---- save -------------------------------------------------------------------
ggsave("accumulation.pdf", p_fig, width = 9.0, height = 3.8)
ggsave("accumulation.png", p_fig, width = 9.0, height = 3.8, dpi = 300)
