## ---- Raw data diagnostics: species overlap across locations (LDB vs RDB) ----
# Reads both MB-LDB and MB-RDB sheets from DryadMetabarcodingData.xlsx,
# extracts observed plant-pollinator pairs per location and method
# using the same logic as prepare_adjacency_matrices.R, and runs
# species overlap diagnostics to investigate community structure.

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(scales)

raw_path  <- "data/raw_data/DryadMetabarcodingData.xlsx"
locations <- c("Zumwalt", "Starkey", "Threemile")

# ---- read and prepare both sheets ----
read_sheet <- function(sheet, db_label) {
  read_excel(raw_path, sheet = sheet) %>%
    mutate(pollinator = str_trim(paste(Genus, Species)),
           database   = db_label)
}

ldb <- read_sheet("MB-LDB", "LDB")
rdb <- read_sheet("MB-RDB", "RDB")
raw <- bind_rows(ldb, rdb)

# ---- extract observed pairs per location × method × database ----
# (same collapse logic as prepare_adjacency_matrices.R: months and sexes pooled)

obs_pairs <- raw %>%
  select(database, location = Location, pollinator, plant = Observation) %>%
  filter(!is.na(plant), plant != "") %>%
  distinct() %>%
  mutate(method = "observation")

mb_pairs <- raw %>%
  select(database, location = Location, pollinator,
         Taxon1, Taxon2, Taxon3, Taxon4, Taxon5,
         Taxon6, Taxon7, Taxon8, Taxon9) %>%
  pivot_longer(cols      = starts_with("Taxon"),
               names_to  = "taxon_col",
               values_to = "plant") %>%
  select(database, location, pollinator, plant) %>%
  filter(!is.na(plant), plant != "") %>%
  distinct() %>%
  mutate(method = "metabarcoding")

df_raw <- bind_rows(obs_pairs, mb_pairs) %>%
  mutate(interaction_id = paste(pollinator, plant, sep = "___"))

# ---- text diagnostics ----
for (db in c("LDB", "RDB")) {
  for (m in c("observation", "metabarcoding")) {
    df_dm <- df_raw %>% filter(database == db, method == m)

    cat("\n==============================\n")
    cat("Database:", db, "| Method:", m, "\n")
    cat("==============================\n")

    # 1. network size per location
    cat("\n--- Network size per location ---\n")
    df_dm %>%
      group_by(location) %>%
      summarise(
        n_plants      = n_distinct(plant),
        n_pollinators = n_distinct(pollinator),
        n_links       = n_distinct(interaction_id),
        .groups = "drop"
      ) %>% print()

    # 2. pairwise species overlap
    cat("\n--- Plant species shared across locations ---\n")
    plants_by_loc <- df_dm %>%
      group_by(location) %>%
      summarise(sp = list(unique(plant)), .groups = "drop")
    locs <- plants_by_loc$location
    for (i in seq_along(locs)) {
      for (j in seq_along(locs)) {
        if (j > i) {
          n <- length(intersect(plants_by_loc$sp[[i]], plants_by_loc$sp[[j]]))
          cat(sprintf("  %s ∩ %s: %d shared plants\n", locs[i], locs[j], n))
        }
      }
    }

    cat("\n--- Pollinator species shared across locations ---\n")
    polls_by_loc <- df_dm %>%
      group_by(location) %>%
      summarise(sp = list(unique(pollinator)), .groups = "drop")
    for (i in seq_along(locs)) {
      for (j in seq_along(locs)) {
        if (j > i) {
          n <- length(intersect(polls_by_loc$sp[[i]], polls_by_loc$sp[[j]]))
          cat(sprintf("  %s ∩ %s: %d shared pollinators\n", locs[i], locs[j], n))
        }
      }
    }

    # 3. link overlap
    cat("\n--- Observed links shared across locations ---\n")
    links_by_loc <- df_dm %>%
      group_by(location) %>%
      summarise(ids = list(unique(interaction_id)), .groups = "drop")
    for (i in seq_along(locs)) {
      for (j in seq_along(locs)) {
        if (j > i) {
          n <- length(intersect(links_by_loc$ids[[i]], links_by_loc$ids[[j]]))
          cat(sprintf("  %s ∩ %s: %d shared links\n", locs[i], locs[j], n))
        }
      }
    }

    # 4. distribution of how many locations each link is observed in
    cat("\n--- Distribution: how many locations each link appears in ---\n")
    df_dm %>%
      group_by(interaction_id) %>%
      summarise(n_locs = n_distinct(location), .groups = "drop") %>%
      count(n_locs) %>%
      print()
  }
}

# ---- community composition tile plot ----
plot_community_raw <- function(df_raw, db_name, method_name) {
  df_dm <- df_raw %>% filter(database == db_name, method == method_name)

  sp_present <- bind_rows(
    df_dm %>% distinct(location, species = plant)       %>% mutate(taxon = "Plants"),
    df_dm %>% distinct(location, species = pollinator)  %>% mutate(taxon = "Pollinators")
  )

  n_locs_per_sp <- sp_present %>%
    group_by(taxon, species) %>%
    summarise(n_locs = n_distinct(location), .groups = "drop")

  plot_df <- sp_present %>%
    distinct(taxon, species) %>%
    tidyr::crossing(location = locations) %>%
    left_join(sp_present %>% mutate(present = TRUE),
              by = c("taxon", "species", "location")) %>%
    left_join(n_locs_per_sp, by = c("taxon", "species")) %>%
    mutate(
      present  = !is.na(present),
      fill_val = factor(ifelse(present, as.character(n_locs), "Absent"),
                        levels = c("Absent", "1", "2", "3"))
    )

  sp_order <- n_locs_per_sp %>%
    arrange(taxon, desc(n_locs), species) %>%
    pull(species)

  plot_df <- plot_df %>%
    mutate(species = factor(species, levels = rev(sp_order)))

  n_species <- n_distinct(plot_df$species)

  ggplot(plot_df, aes(x = location, y = species, fill = fill_val)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_manual(
      values = c("Absent" = "grey92", "1" = "#aec6cf", "2" = "#f4a261", "3" = "#2a9d8f"),
      labels = c("Absent", "1 location", "2 locations", "3 locations"),
      name   = "Presence"
    ) +
    facet_wrap(~ taxon, scales = "free_y", ncol = 2) +
    labs(title   = paste("Community composition —", db_name, "|", method_name),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y     = element_text(size = max(4, min(8, 400 / n_species)), face = "italic"),
      axis.text.x     = element_text(size = 11, face = "bold"),
      strip.text      = element_text(size = 12, face = "bold"),
      panel.grid      = element_blank(),
      legend.position = "bottom",
      plot.title      = element_text(face = "bold")
    )
}

for (db in c("LDB", "RDB")) {
  for (m in c("observation", "metabarcoding")) {
    print(plot_community_raw(df_raw, db, m))
  }
}

# ---- interaction richness column plots (locations side by side) ----
make_richness_plot_raw <- function(df_raw, db_name, method_name,
                                   location_name, taxon_name) {
  df_dm <- df_raw %>%
    filter(database == db_name, method == method_name, location == location_name)

  richness <- if (taxon_name == "Plants") {
    df_dm %>%
      group_by(species = plant) %>%
      summarise(n_partners = n_distinct(pollinator), .groups = "drop")
  } else {
    df_dm %>%
      group_by(species = pollinator) %>%
      summarise(n_partners = n_distinct(plant), .groups = "drop")
  }

  richness <- richness %>%
    mutate(species = factor(species, levels = sort(unique(species), decreasing = TRUE)))

  ggplot(richness, aes(x = n_partners, y = species)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n_partners), hjust = -0.3, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title    = paste(location_name, "—", taxon_name),
         subtitle = paste(db_name, "|", method_name),
         x        = "Number of interaction partners",
         y        = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y        = element_text(size = 8, face = "italic"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(face = "bold"),
      plot.subtitle      = element_text(color = "grey40")
    )
}

for (db in c("LDB", "RDB")) {
  for (m in c("observation", "metabarcoding")) {
    for (tx in c("Plants", "Pollinators")) {
      combined <- make_richness_plot_raw(df_raw, db, m, locations[1], tx) |
                  make_richness_plot_raw(df_raw, db, m, locations[2], tx) |
                  make_richness_plot_raw(df_raw, db, m, locations[3], tx)
      print(combined +
              plot_annotation(title = paste(tx, "—", db, "|", m),
                              theme = theme(plot.title = element_text(face = "bold",
                                                                      size = 13))))
    }
  }
}

# not enough overlap between communities!! :(((