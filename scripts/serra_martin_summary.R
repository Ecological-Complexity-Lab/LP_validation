library(tidyverse)
library(vegan)

# ── 1. Load & clean ──────────────────────────────────────────────────────────
df_raw <- read_delim(
  "data/raw_data/serra_martin_pollination/cabrera_22_23_habitat.csv",
  delim = ";",
  show_col_types = FALSE
)

df <- df_raw %>%
  # Harmonise plant names (following the authors' Python scripts)
  mutate(`Plant sp` = recode(`Plant sp`,
    "Daucus carota L. subsp. Majoricus" = "Daucus carota",
    "Rosmarinus officinalis"            = "Salvia rosmarinus"
  )) %>%
  # Drop non-standard plant species (orchids used for camera training only)
  filter(!`Plant sp` %in% c("Anacamptis pyramidalis", "Gladiolus communis")) %>%
  # Drop zero-interaction census rows (Pollinator is empty)
  filter(!is.na(Pollinator), Pollinator != "") %>%
  # Derive helper columns
  mutate(
    habitat_type = str_remove(habitat, " [12]$"),   # broad type (3 levels)
    interaction  = paste(`Plant sp`, Pollinator, sep = " × ")
  )

cat("=== DATASET OVERVIEW ===\n")
cat(sprintf("Total interaction records : %d\n", nrow(df)))
cat(sprintf("Sampling years            : %s\n",
            paste(sort(unique(year(dmy_hm(df$Start_T)))), collapse = ", ")))
cat(sprintf("Visits (rounds)           : %d\n", n_distinct(df$visita)))
cat(sprintf("Censuses                  : %d\n", n_distinct(df$censo)))
cat(sprintf("Methods                   : %s\n",
            paste(sort(unique(df$Method)), collapse = ", ")))
cat(sprintf("Habitat patches           : %s\n",
            paste(sort(unique(df$habitat)), collapse = ", ")))
cat("\n")


# ── 2. Helper: richness summary ───────────────────────────────────────────────
richness_summary <- function(data, group_var) {
  data %>%
    group_by({{ group_var }}) %>%
    summarise(
      n_censuses       = n_distinct(censo),
      n_visits         = n_distinct(visita),
      n_plant_spp      = n_distinct(`Plant sp`),
      n_pollinator_spp = n_distinct(Pollinator),
      n_interactions   = n_distinct(interaction),
      total_ind        = sum(`N ind`, na.rm = TRUE),
      .groups = "drop"
    )
}


# ── 3. Summary by habitat PATCH (6 sites) ────────────────────────────────────
cat("=== SUMMARY BY HABITAT PATCH ===\n")
patch_summary <- richness_summary(df, habitat)
print(patch_summary, n = Inf)
cat("\n")


# ── 4. Summary by habitat TYPE (3 broad categories) ──────────────────────────
cat("=== SUMMARY BY HABITAT TYPE ===\n")
type_summary <- richness_summary(df, habitat_type)
print(type_summary, n = Inf)
cat("\n")


# ── 5. Summary by sampling METHOD ────────────────────────────────────────────
cat("=== SUMMARY BY SAMPLING METHOD ===\n")
method_summary <- richness_summary(df, Method) %>%
  mutate(Method = recode(Method,
    "obs" = "Direct observation",
    "rpi" = "Raspberry Pi camera (ACS)"
  ))
print(method_summary, n = Inf)
cat("\n")

# Method × habitat type cross-tab
cat("--- Method × habitat type ---\n")
method_habitat <- df %>%
  group_by(habitat_type, Method) %>%
  summarise(
    n_plant_spp      = n_distinct(`Plant sp`),
    n_pollinator_spp = n_distinct(Pollinator),
    n_interactions   = n_distinct(interaction),
    .groups = "drop"
  ) %>%
  mutate(Method = recode(Method,
    "obs" = "obs (direct)",
    "rpi" = "rpi (camera)"
  )) %>%
  pivot_wider(
    names_from  = Method,
    values_from = c(n_plant_spp, n_pollinator_spp, n_interactions)
  )
print(method_habitat, n = Inf)
cat("\n")


# ── 6. Pairwise overlap between habitat PATCHES ───────────────────────────────
# Using Jaccard similarity (presence/absence) and Bray-Curtis dissimilarity
# (abundance-weighted), computed separately for plant species, pollinator
# species, and interactions.

make_community_matrix <- function(data, group_var, species_var,
                                  abund_col = NULL) {
  g_name <- as_label(enquo(group_var))
  if (is.null(abund_col)) {
    data %>%
      distinct({{ group_var }}, {{ species_var }}) %>%
      mutate(presence = 1L) %>%
      pivot_wider(names_from = {{ species_var }},
                  values_from = presence, values_fill = 0L) %>%
      column_to_rownames(g_name)
  } else {
    data %>%
      group_by({{ group_var }}, {{ species_var }}) %>%
      summarise(abund = sum(.data[[abund_col]], na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = {{ species_var }},
                  values_from = abund, values_fill = 0) %>%
      column_to_rownames(g_name)
  }
}

pairwise_jaccard <- function(mat) {
  d  <- vegdist(mat, method = "jaccard", binary = TRUE)
  sim <- 1 - as.matrix(d)
  sim[lower.tri(sim, diag = TRUE)] <- NA
  as.data.frame(sim) %>%
    rownames_to_column("patch_A") %>%
    pivot_longer(-patch_A, names_to = "patch_B", values_to = "jaccard_sim") %>%
    filter(!is.na(jaccard_sim)) %>%
    arrange(desc(jaccard_sim))
}

pairwise_braycurtis <- function(mat) {
  d <- vegdist(mat, method = "bray")
  sim <- 1 - as.matrix(d)
  sim[lower.tri(sim, diag = TRUE)] <- NA
  as.data.frame(sim) %>%
    rownames_to_column("patch_A") %>%
    pivot_longer(-patch_A, names_to = "patch_B", values_to = "bray_sim") %>%
    filter(!is.na(bray_sim)) %>%
    arrange(desc(bray_sim))
}

# Build matrices by habitat PATCH
plant_mat_patch  <- make_community_matrix(df, habitat, `Plant sp`)
poll_mat_patch   <- make_community_matrix(df, habitat, Pollinator)
inter_mat_patch  <- make_community_matrix(df, habitat, interaction)

plant_abund_mat_patch <- make_community_matrix(df, habitat, `Plant sp`,
                                               abund_col = "N open flowers")
poll_abund_mat_patch  <- make_community_matrix(df, habitat, Pollinator,
                                               abund_col = "N ind")

cat("=== PAIRWISE OVERLAP BETWEEN HABITAT PATCHES ===\n")
cat("\n-- Plant species Jaccard similarity (presence/absence) --\n")
print(pairwise_jaccard(plant_mat_patch))

cat("\n-- Pollinator species Jaccard similarity (presence/absence) --\n")
print(pairwise_jaccard(poll_mat_patch))

cat("\n-- Interaction Jaccard similarity (presence/absence) --\n")
print(pairwise_jaccard(inter_mat_patch))

cat("\n-- Plant species Bray-Curtis similarity (abundance-weighted) --\n")
print(pairwise_braycurtis(plant_abund_mat_patch))

cat("\n-- Pollinator species Bray-Curtis similarity (abundance-weighted) --\n")
print(pairwise_braycurtis(poll_abund_mat_patch))


# ── 7. Pairwise overlap between habitat TYPES ─────────────────────────────────
plant_mat_type  <- make_community_matrix(df, habitat_type, `Plant sp`)
poll_mat_type   <- make_community_matrix(df, habitat_type, Pollinator)
inter_mat_type  <- make_community_matrix(df, habitat_type, interaction)

plant_abund_mat_type <- make_community_matrix(df, habitat_type, `Plant sp`,
                                              abund_col = "N open flowers")
poll_abund_mat_type  <- make_community_matrix(df, habitat_type, Pollinator,
                                              abund_col = "N ind")

cat("\n=== PAIRWISE OVERLAP BETWEEN HABITAT TYPES ===\n")

cat("\n-- Plant species Jaccard similarity --\n")
print(pairwise_jaccard(plant_mat_type))

cat("\n-- Pollinator species Jaccard similarity --\n")
print(pairwise_jaccard(poll_mat_type))

cat("\n-- Interaction Jaccard similarity --\n")
print(pairwise_jaccard(inter_mat_type))

cat("\n-- Plant species Bray-Curtis similarity (abundance-weighted) --\n")
print(pairwise_braycurtis(plant_abund_mat_type))

cat("\n-- Pollinator species Bray-Curtis similarity (abundance-weighted) --\n")
print(pairwise_braycurtis(poll_abund_mat_type))


# ── 8. Species and interaction lists per habitat type ─────────────────────────
# (useful as a quick reference for what is exclusive to each habitat)
cat("\n=== EXCLUSIVE SPECIES PER HABITAT TYPE ===\n")

exclusive_species <- function(data, group_var, species_var) {
  g  <- enquo(group_var)
  sp <- enquo(species_var)
  sp_by_group <- data %>%
    group_by(!!g) %>%
    summarise(species = list(unique(!!sp)), .groups = "drop")

  map_dfr(seq_len(nrow(sp_by_group)), function(i) {
    focal_group   <- sp_by_group[[as_label(g)]][i]
    focal_species <- sp_by_group$species[[i]]
    other_species <- unlist(sp_by_group$species[-i])
    exclusive     <- setdiff(focal_species, other_species)
    tibble(
      group     = focal_group,
      n_total   = length(focal_species),
      n_exclusive = length(exclusive),
      pct_exclusive = round(100 * length(exclusive) / length(focal_species), 1),
      exclusive_spp = paste(sort(exclusive), collapse = "; ")
    )
  })
}

cat("\n-- Exclusive PLANT species per habitat type --\n")
plant_excl <- exclusive_species(df, habitat_type, `Plant sp`)
print(plant_excl %>% select(-exclusive_spp), n = Inf)
for (i in seq_len(nrow(plant_excl))) {
  cat(sprintf("\n  %s exclusive plants: %s\n",
              plant_excl$group[i], plant_excl$exclusive_spp[i]))
}

cat("\n-- Exclusive POLLINATOR species per habitat type --\n")
poll_excl <- exclusive_species(df, habitat_type, Pollinator)
print(poll_excl %>% select(-exclusive_spp), n = Inf)
for (i in seq_len(nrow(poll_excl))) {
  cat(sprintf("\n  %s exclusive pollinators: %s\n",
              poll_excl$group[i], poll_excl$exclusive_spp[i]))
}

cat("\n-- Exclusive INTERACTIONS per habitat type --\n")
inter_excl <- exclusive_species(df, habitat_type, interaction)
print(inter_excl %>% select(-exclusive_spp), n = Inf)
