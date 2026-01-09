library(dplyr)
library(eulerr)
library(ggplot2)

# ---- observed links ----
## ---- euler diagram ----
# 1) De-duplicate (one row per interaction_id × method)
df0 <- df_categorized %>%
  mutate(method = recode(method,
                         "method1" = "mist_nets",
                         "method2" = "observations")) %>%
  group_by(method, interaction_id) %>%
  slice(1) %>%                 # keep one record per method×interaction
  ungroup()

methods <- df0 %>% distinct(method) %>% arrange(method) %>% pull(method)
stopifnot(length(methods) == 2)
m1 <- methods[1]
m2 <- methods[2]

# 2) Observed interaction sets per method
A <- df0 %>% filter(method == m1, original_binary == 1) %>% pull(interaction_id) %>% unique()
B <- df0 %>% filter(method == m2, original_binary == 1) %>% pull(interaction_id) %>% unique()

A_only <- setdiff(A, B)
B_only <- setdiff(B, A)
AB     <- intersect(A, B)

# 3) Fit Euler from sets (robust)
fit <- eulerr::euler(setNames(list(A, B), c(m1, m2)))

# 4) Category breakdown counts for labels
get_n <- function(df, cat) {
  out <- df %>% filter(link_category == cat) %>% summarise(n = sum(n)) %>% pull(n)
  ifelse(length(out) == 0 || is.na(out), 0, out)
}

m1_only_counts <- df0 %>%
  filter(method == m1, interaction_id %in% A_only,
         link_category %in% c("locally_unique_links", "unsupported_links")) %>%
  count(link_category, name = "n")

m2_only_counts <- df0 %>%
  filter(method == m2, interaction_id %in% B_only,
         link_category %in% c("locally_unique_links", "unsupported_links")) %>%
  count(link_category, name = "n")

# overlap categories appear once per method -> divide by 2
ab_counts <- df0 %>%
  filter(interaction_id %in% AB,
         link_category %in% c("confirmed_links", "cryptic_links")) %>% # cryptic links should not be here, they are counted for each method seperately
  count(link_category, name = "n") %>% # for recurrent links we need to count the unique interactions.
  mutate(n = n / 2)

lab_A <- paste0(
  "Method-specific observed\n",
  "locally_unique: ", get_n(m1_only_counts, "locally_unique_links"), "\n",
  "unsupported: ",   get_n(m1_only_counts, "unsupported_links")
)

lab_B <- paste0(
  "Method-specific observed\n",
  "locally_unique: ", get_n(m2_only_counts, "locally_unique_links"), "\n",
  "unsupported: ",   get_n(m2_only_counts, "unsupported_links")
)

lab_AB <- paste0(
  "Observed in both\n",
  "confirmed: ", get_n(ab_counts, "confirmed_links"), "\n",
  "cryptic: ",   get_n(ab_counts, "cryptic_links")
)

# 5) Get ellipse centers for placing labels
# fit$ellipses is a data.frame-like object with x,y as center coordinates
ell <- as.data.frame(fit$ellipses)
# rownames are the set names
cx1 <- ell[m1, "x"]; cy1 <- ell[m1, "y"]
cx2 <- ell[m2, "x"]; cy2 <- ell[m2, "y"]

# Put overlap label roughly between centers (works well for 2-set Venn)
cx_ab <- (cx1 + cx2) / 2
cy_ab <- (cy1 + cy2) / 2

# 6) Plot and annotate
p <- plot(
  fit,
  fills = list(fill = c("#4C78A8", "#F58518"), alpha = 0.35),
  edges = list(col = "white", lwd = 1),
  labels = list(font = 2)
)

p

## ---- ggplot with count labels ----

# ---- 1) Observed sets per method (unique interaction IDs) ---

# helper: count unique interactions in a subset + category
n_uid <- function(df, cat) df %>% filter(link_category == cat) %>% summarise(n = n_distinct(interaction_id)) %>% pull(n)

# ---- 2) Method-only observed category counts (unique IDs) ---
m1_only_df <- df0 %>%
  filter(method == m1, interaction_id %in% A_only, original_binary == 1)

m2_only_df <- df0 %>%
  filter(method == m2, interaction_id %in% B_only, original_binary == 1)

m1_loc <- n_uid(m1_only_df, "locally_unique_links")
m1_uns <- n_uid(m1_only_df, "unsupported_links")

m2_loc <- n_uid(m2_only_df, "locally_unique_links")
m2_uns <- n_uid(m2_only_df, "unsupported_links")

# ---- 3) Overlap (shared observed) counts + cryptic split (unique IDs) ---
# Build one row per interaction_id with both methods' category/prediction
wide_AB <- df0 %>%
  filter(interaction_id %in% AB) %>%
  select(interaction_id, method, predicted_bin, link_category) %>%
  pivot_wider(
    names_from = method,
    values_from = c(predicted_bin, link_category)
  )

# totals in overlap (unique interactions)
confirmed_total <- sum(
  wide_AB[[paste0("link_category_", m1)]] == "confirmed_links" &
    wide_AB[[paste0("link_category_", m2)]] == "confirmed_links",
  na.rm = TRUE
)

cryptic_total <- sum(
  wide_AB[[paste0("link_category_", m1)]] == "cryptic_links" &
    wide_AB[[paste0("link_category_", m2)]] == "cryptic_links",
  na.rm = TRUE
)

# NOTE: if an interaction is shared-observed, its per-method category should be either confirmed or cryptic (by your rules).
# But to be safe, we compute cryptic split from per-method categories:
cryptic_only_m1 <- sum(
  wide_AB[[paste0("link_category_", m1)]] == "cryptic_links" &
    wide_AB[[paste0("link_category_", m2)]] != "cryptic_links",
  na.rm = TRUE
)

cryptic_only_m2 <- sum(
  wide_AB[[paste0("link_category_", m2)]] == "cryptic_links" &
    wide_AB[[paste0("link_category_", m1)]] != "cryptic_links",
  na.rm = TRUE
)

cryptic_both <- sum(
  wide_AB[[paste0("link_category_", m1)]] == "cryptic_links" &
    wide_AB[[paste0("link_category_", m2)]] == "cryptic_links",
  na.rm = TRUE
)

# You may also want confirmed split (usually all confirmed should be confirmed in both)
# confirmed_total is computed as both-confirmed by design.

# ---- 4) Labels (all are unique ID counts) ---
lab_A <- paste0(
  "Method-specific observed\n",
  "locally unique: ", m1_loc, "\n",
  "unsupported: ", m1_uns, "\n",
  "Total: ", length(A_only)
)

lab_B <- paste0(
  "Method-specific observed\n",
  "locally unique: ", m2_loc, "\n",
  "unsupported: ", m2_uns, "\n",
  "Total: ", length(B_only)
)

lab_AB <- paste0(
  "Observed in both\n",
  "confirmed: ", confirmed_total, "\n",
  "cryptic: ", (cryptic_only_m1 + cryptic_only_m2 + cryptic_both), "\n\n",
  "Cryptic split:\n",
  m1, " only: ", cryptic_only_m1, "\n",
  m2, " only: ", cryptic_only_m2, "\n",
  "both: ", cryptic_both, "\n",
  "Total: ", length(AB)
)

# ---- 5) Draw circles manually (ggplot2) ---
circle_df <- function(cx, cy, r, n = 400) {
  t <- seq(0, 2*pi, length.out = n)
  data.frame(x = cx + r*cos(t), y = cy + r*sin(t))
}

r <- 1.8
c1 <- c(-1.2, 0)   # center for m1
c2 <- c( 1.2, 0)   # center for m2

circ1 <- circle_df(c1[1], c1[2], r) %>% mutate(set = m1)
circ2 <- circle_df(c2[1], c2[2], r) %>% mutate(set = m2)

pos <- data.frame(
  x = c(c1[1] - 0.6, 0, c2[1] + 0.6),
  y = c(0, 0, 0),
  label = c(lab_A, lab_AB, lab_B)
)

names_pos <- data.frame(
  x = c(c1[1], c2[1]),
  y = c(r + 0.35, r + 0.35),
  label = c(m1, m2)
)

ggplot() +
  geom_polygon(data = circ1, aes(x, y), fill = "#4C78A8", alpha = 0.35, color = "white", linewidth = 0.8) +
  geom_polygon(data = circ2, aes(x, y), fill = "#F58518", alpha = 0.35, color = "white", linewidth = 0.8) +
  geom_text(data = names_pos, aes(x, y, label = label), fontface = "bold", size = 4) +
  geom_text(data = pos, aes(x, y, label = label), size = 3.4, lineheight = 1.05) +
  coord_equal() +
  theme_void(base_size = 12) +
  labs(title = "Observed interactions: overlap between methods") +
  theme(plot.title = element_text(face = "bold"))

## ---- rectangles ----
# --- counts (unique interaction IDs) ---
A1  <- length(A)      # observed in mist_nets
A2  <- length(B)      # observed in observations
A12 <- length(AB)     # observed in both

stopifnot(A12 <= A1, A12 <= A2)

# --- scale factor: controls overall size of the figure (bigger = larger rectangles) ---
s <- 0.5  # try 0.5 or 2 if you want smaller/larger

# --- choose square-like rectangles so area = count*s ---
w1 <- sqrt(A1  * s); h1 <- sqrt(A1  * s)
w2 <- sqrt(A2  * s); h2 <- sqrt(A2  * s)
wo <- sqrt(A12 * s); ho <- sqrt(A12 * s)

# ensure overlap fits inside both rectangles
stopifnot(wo <= w1, ho <= h1, wo <= w2, ho <= h2)

# --- rectangle coordinates ---
# Rect 1 anchored at origin
r1 <- data.frame(
  method = m1,
  xmin = 0, xmax = w1,
  ymin = 0, ymax = h1
)

# Rect 2 placed so its bottom-left corner overlaps rect1's top-right corner
r2 <- data.frame(
  method = m2,
  xmin = w1 - wo, xmax = (w1 - wo) + w2,
  ymin = h1 - ho, ymax = (h1 - ho) + h2
)

rects <- rbind(r1, r2)

# --- label positions (tweak these multipliers if needed) ---
lab_pos <- data.frame(
  region = c("A_only", "AB", "B_only"),
  x = c(w1 * 0.30, w1 - wo/2, (w1 - wo) + w2 * 0.70),
  y = c(h1 * 0.30, h1 - ho/2, (h1 - ho) + h2 * 0.70),
  label = c(lab_A, lab_AB, lab_B)
)

method_pos <- data.frame(
  x = c(w1 * 0.20, (w1 - wo) + w2 * 0.80),
  y = c(h1 + 0.10*h1, (h1 - ho) + h2 + 0.10*h2),
  label = c(m1, m2)
)

# --- plot ---
ggplot() +
  geom_rect(data = rects,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = method),
            alpha = 0.30, color = "white", linewidth = 0.9) +
  geom_text(data = method_pos, aes(x, y, label = label),
            fontface = "bold", size = 4) +
  geom_text(data = lab_pos, aes(x, y, label = label),
            size = 3.4, lineheight = 1.05) +
  coord_equal() +
  theme_void(base_size = 12) +
  labs(title = "Observed interactions: overlap between methods (area ∝ count)") +
  theme(plot.title = element_text(face = "bold")) +
  scale_fill_manual(values = c("#4C78A8", "#F58518"), guide = "none")

# ---- unobserved links ----
# ---- 1) Define UNOBSERVED sets per method (unique interaction IDs) ---
U1 <- df0 %>% filter(method == m1, original_binary == 0) %>% pull(interaction_id) %>% unique()
U2 <- df0 %>% filter(method == m2, original_binary == 0) %>% pull(interaction_id) %>% unique()

U1_only <- setdiff(U1, U2)     # unobserved only in m1 (observed in m2)
U2_only <- setdiff(U2, U1)     # unobserved only in m2 (observed in m1)
U12     <- intersect(U1, U2)   # unobserved in both

# ---- 2) Count categories by unique interaction IDs ---
n_uid <- function(df, cat) df %>% filter(link_category == cat) %>%
  summarise(n = n_distinct(interaction_id)) %>% pull(n)

# Non-overlap (method-specific unobserved regions)
m1_unobs_only_df <- df0 %>%
  filter(method == m1, interaction_id %in% U1_only, original_binary == 0)

m2_unobs_only_df <- df0 %>%
  filter(method == m2, interaction_id %in% U2_only, original_binary == 0)

m1_feas <- n_uid(m1_unobs_only_df, "feasible_links")
m1_miss <- n_uid(m1_unobs_only_df, "possibly_missing_links")

m2_feas <- n_uid(m2_unobs_only_df, "feasible_links")
m2_miss <- n_uid(m2_unobs_only_df, "possibly_missing_links")

# Overlap: unobserved in both (compute from wide table so it's unique IDs)
wide_U12 <- df0 %>%
  filter(interaction_id %in% U12) %>%
  select(interaction_id, method, link_category, predicted_bin) %>%
  pivot_wider(names_from = method, values_from = c(link_category, predicted_bin))

# Totals in overlap (unique interactions)
forbidden_total <- sum(
  wide_U12[[paste0("link_category_", m1)]] == "likely_forbidden" &
    wide_U12[[paste0("link_category_", m2)]] == "likely_forbidden",
  na.rm = TRUE
)

spurious_total <- sum(
  wide_U12[[paste0("link_category_", m1)]] == "spurious_links" &
    wide_U12[[paste0("link_category_", m2)]] == "spurious_links",
  na.rm = TRUE
)

# (Optional) split spurious by method if you want symmetry with the cryptic split:
spurious_only_m1 <- sum(
  wide_U12[[paste0("link_category_", m1)]] == "spurious_links" &
    wide_U12[[paste0("link_category_", m2)]] != "spurious_links",
  na.rm = TRUE
)
spurious_only_m2 <- sum(
  wide_U12[[paste0("link_category_", m2)]] == "spurious_links" &
    wide_U12[[paste0("link_category_", m1)]] != "spurious_links",
  na.rm = TRUE
)
spurious_both <- sum(
  wide_U12[[paste0("link_category_", m1)]] == "spurious_links" &
    wide_U12[[paste0("link_category_", m2)]] == "spurious_links",
  na.rm = TRUE
)

# ---- 3) Build labels (all unique ID counts) ---
lab_U1 <- paste0(
  "Unobserved here,\nobserved elsewhere\n",
  "feasible: ", m1_feas, "\n",
  "missing: ",  m1_miss, "\n",
  "Total: ", length(U1_only)
)

lab_U2 <- paste0(
  "Unobserved here,\nobserved elsewhere\n",
  "feasible: ", m2_feas, "\n",
  "missing: ",  m2_miss, "\n",
  "Total: ", length(U2_only)
)

lab_U12 <- paste0(
  "Unobserved in both\n",
  "likely forbidden: ", forbidden_total, "\n",
  "spurious: ", (spurious_only_m1 + spurious_only_m2 + spurious_both), "\n\n",
  "Spurious split:\n",
  m1, " only: ", spurious_only_m1, "\n",
  m2, " only: ", spurious_only_m2, "\n",
  "both: ", spurious_both, "\n",
  "Total: ", length(U12)
)

# ---- 4) Rectangles with corner overlap (area ∝ count) ----
U1n  <- length(U1)
U2n  <- length(U2)
U12n <- length(U12)

stopifnot(U12n <= U1n, U12n <= U2n)

s <- 0.5  # scale factor (change if you want overall size)

w1 <- sqrt(U1n  * s); h1 <- sqrt(U1n  * s)
w2 <- sqrt(U2n  * s); h2 <- sqrt(U2n  * s)
wo <- sqrt(U12n * s); ho <- sqrt(U12n * s)

stopifnot(wo <= w1, ho <= h1, wo <= w2, ho <= h2)

r1 <- data.frame(method = m1, xmin = 0, xmax = w1, ymin = 0, ymax = h1)
r2 <- data.frame(method = m2,
                 xmin = w1 - wo, xmax = (w1 - wo) + w2,
                 ymin = h1 - ho, ymax = (h1 - ho) + h2)

rects <- rbind(r1, r2)

lab_pos <- data.frame(
  x = c(w1 * 0.30, w1 - wo/2, (w1 - wo) + w2 * 0.70),
  y = c(h1 * 0.30, h1 - ho/2, (h1 - ho) + h2 * 0.70),
  label = c(lab_U1, lab_U12, lab_U2)
)

method_pos <- data.frame(
  x = c(w1 * 0.20, (w1 - wo) + w2 * 0.80),
  y = c(h1 + 0.10*h1, (h1 - ho) + h2 + 0.10*h2),
  label = c(m1, m2)
)

ggplot() +
  geom_rect(
    data = rects,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = method),
    alpha = 0.30, color = "white", linewidth = 0.9
  ) +
  geom_text(data = method_pos, aes(x, y, label = label),
            fontface = "bold", size = 4) +
  geom_text(data = lab_pos, aes(x, y, label = label),
            size = 3.4, lineheight = 1.05) +
  coord_equal() +
  theme_void(base_size = 12) +
  labs(title = "Unobserved interactions: overlap between methods (area ∝ count)") +
  theme(plot.title = element_text(face = "bold")) +
  scale_fill_manual(values = c("thistle3", "peachpuff3"), guide = "none")
