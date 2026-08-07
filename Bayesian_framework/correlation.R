## ---------------------------------------------------------------------------
## SI figure: correlated errors between the local observation and the replicates
## Evidence throughout: Yhat = 1, O_l = 0, detection in every replicate (n = R)
## ---------------------------------------------------------------------------

library(ggplot2)
library(patchwork)

eY     <- 0.20      # model error
rho    <- 0.15      # realisation rate
el_fix <- 0.30      # sampling miss rate, when held fixed
f_fix  <- 0.05      # false-detection rate, when held fixed
nu     <- 40        # concentration of the beta distribution over the shared rate
Rs     <- 1:50      # number of replicates

## the eight categories, as (z_Y, z_l, z_r) signatures
sig <- rbind(c(1,1,1), c(1,1,0), c(1,0,1), c(1,0,0),
             c(0,1,1), c(0,1,0), c(0,0,1), c(0,0,0))
rownames(sig) <- c("recurrent", "locally unique", "possibly missing", "phantom",
                   "model-elusive", "weakly-supported", "locally absent",
                   "possibly forbidden")
pm <- 3   # row of possibly missing

## Likelihood of the evidence for each category, vectorised over a grid of rates.
## Directional local axis throughout (Eq. likelihood_dir): a realised link is
## missed at epsilon_l, an absent one is recorded at f.
lik <- function(el, f, R) {
  n  <- max(length(el), length(f))
  el <- rep_len(el, n); f <- rep_len(f, n)
  lmatch <- 1 - f                                   # O_l = 0 matches z_l = 0
  p1 <- rho * (1 - el)
  out <- matrix(0, nrow = 8, ncol = n)
  for (i in 1:8) {
    lY <- if (sig[i, 1] == 1) 1 - eY else eY
    lL <- if (sig[i, 2] == 0) lmatch else el
    p  <- if (sig[i, 3] == 1) p1 else f
    out[i, ] <- lY * lL * p^R
  }
  out
}

## posterior with both rates fixed
post_fixed <- function(el, f, R) {
  L <- lik(el, f, R); as.vector(L / sum(L))
}

## posterior after marginalising over a rate drawn once per link
grid <- (seq_len(20000) - 0.5) / 20000
post_shared <- function(R, w, el = el_fix, f = f_fix) {
  tot <- as.vector(lik(el, f, R) %*% w)
  tot / sum(tot)
}

## ---- (a) detectability u = 1 - epsilon_l drawn once per link ----------------
u_bar <- 1 - el_fix
al    <- u_bar * nu ; bl <- (1 - u_bar) * nu
w_u   <- dbeta(grid, al, bl)
a_ind <- sapply(Rs, function(R) post_fixed(el_fix, f_fix, R)[pm])
a_sha <- sapply(Rs, function(R) post_shared(R, w_u, el = 1 - grid, f = f_fix)[pm])

## closed form of the two replacements in Eq. (corr_el), valid when n = R
post_closed <- function(n) {
  prod_u <- prod((al + seq_len(n) - 1) / (al + bl + seq_len(n) - 1))
  out <- numeric(8)
  for (i in 1:8) {
    lY   <- if (sig[i, 1] == 1) 1 - eY else eY
    el_t <- bl / (al + bl + n * sig[i, 3])              # miss rate, updated only if z_r = 1
    lL   <- if (sig[i, 2] == 0) 1 - f_fix else el_t     # absent locally: f, untouched
    Pc   <- if (sig[i, 3] == 1) rho^n * prod_u else f_fix^n
    out[i] <- lY * lL * Pc
  }
  out / sum(out)
}
stopifnot(max(abs(a_sha - sapply(Rs, function(R) post_closed(R)[pm]))) < 1e-6)

## ---- (b) false-detection rate f drawn once per link -------------------------
w_f   <- dbeta(grid, f_fix * nu, (1 - f_fix) * nu)
b_ind <- sapply(Rs, function(R) post_fixed(el_fix, f_fix, R)[pm])
b_sha <- sapply(Rs, function(R) post_shared(R, w_f, el = el_fix, f = grid)[pm])

kappa <- (1 - eY) * (1 - f_fix) / ((1 - f_fix) + el_fix)   # 0.608

dat <- rbind(
  data.frame(panel = "a", R = Rs, errors = "independent", p = a_ind),
  data.frame(panel = "a", R = Rs, errors = "shared",      p = a_sha),
  data.frame(panel = "b", R = Rs, errors = "independent", p = b_ind),
  data.frame(panel = "b", R = Rs, errors = "shared",      p = b_sha))
dat$errors <- factor(dat$errors, levels = c("independent", "shared"))

cols <- c(independent = "#2a78d6", shared = "#eb6834")
pct  <- function(x) paste0(round(100 * x), "%")

base <- theme_classic(base_size = 10) +
  theme(legend.position   = "bottom",
        legend.title      = element_blank(),
        legend.key.width  = unit(18, "pt"),
        plot.title        = element_text(size = 10, hjust = 0),
        axis.line         = element_line(colour = "grey40", linewidth = 0.3),
        axis.ticks        = element_line(colour = "grey40", linewidth = 0.3))

## the dashed rule is drawn AFTER the curves: the independent curve rises to
## kappa, so a rule drawn underneath it is hidden by the blue line at large R
pa <- ggplot(subset(dat, panel == "a"), aes(R, p, colour = errors)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = kappa, linetype = "22",
             colour = "grey25", linewidth = 0.45) +
  annotate("text", x = max(Rs), y = kappa, label = "kappa == 0.608", parse = TRUE,
           hjust = 1, vjust = 1.9, size = 3, colour = "grey25") +
  scale_colour_manual(values = cols) +
  scale_y_continuous(limits = c(0.40, 0.82), labels = pct) +
  labs(title = "(a) sampling miss rate shared",
       x = "Replicates, R", y = "P(possibly missing | E)") + base

pb <- ggplot(subset(dat, panel == "b"), aes(R, p, colour = errors)) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = cols) +
  scale_y_continuous(limits = c(0, 0.66), labels = pct) +
  labs(title = "(b) false-detection rate shared",
       x = "Replicates, R", y = NULL) + base

fig <- pa + pb + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("correlation.pdf", fig, width = 7.2, height = 3.2)
