## ---- Raw data diagnostics ----
# Section 1: Species overlap across locations (LDB vs RDB)
#   Reads both MB-LDB and MB-RDB sheets from DryadMetabarcodingData.xlsx,
#   extracts observed plant-pollinator pairs per location and method
#   using the same logic as prepare_adjacency_matrices.R, and runs
#   species overlap diagnostics to investigate community structure.
#
# Section 2: EUPPollNet multi-location multi-method study search
#   Reads the EUPPollNet database (Interaction_data.rds + Metadata.csv) and
#   identifies studies that (a) sampled the same site across multiple visits or
#   (b) belong to a cluster of studies sharing a geographic location — with at
#   least two different sampling methods represented in that cluster.

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
# the sites are within different habitats.

# =============================================================================
# SECTION 2: EUPPollNet — find studies with multiple sites AND multiple methods
# Goal: identify Study_id entries that have ≥ 2 distinct Network_id values
# AND ≥ 2 distinct Sampling_method values (from Metadata).
# =============================================================================

eup_int_path  <- "data/raw_data/euppollnet/Interaction_data.rds"
eup_meta_path <- "data/raw_data/euppollnet/Metadata.csv"

eup_int  <- readRDS(eup_int_path)
eup_meta <- read.csv(eup_meta_path, stringsAsFactors = FALSE)

# ---- peek at structure ----
cat("\n==============================\n")
cat("EUPPollNet: data structure\n")
cat("==============================\n")
cat("Interaction data — dimensions:", paste(dim(eup_int), collapse = " x "),
    "| columns:\n")
print(names(eup_int))
cat("\nFirst 3 rows:\n")
print(head(eup_int, 3))
cat("\nMetadata — dimensions:", paste(dim(eup_meta), collapse = " x "),
    "| columns:\n")
print(names(eup_meta))
cat("\nUnique Sampling_method values:\n")
print(sort(unique(eup_meta$Sampling_method)))

# ---- attach Sampling_method if not already in the interaction data ----
if ("Sampling_method" %in% names(eup_int)) {
  eup <- eup_int
} else {
  eup <- eup_int %>%
    left_join(eup_meta %>% select(Study_id, Sampling_method), by = "Study_id")
}

# ---- one row per Network_id with coordinates and study metadata ----
# Adjust lat_col / lon_col if the actual column names differ (see peek output).
lat_col <- "Latitude"
lon_col <- "Longitude"

if (!lat_col %in% names(eup) || !lon_col %in% names(eup)) {
  stop("Coordinate columns '", lat_col, "' / '", lon_col,
       "' not found. Update lat_col / lon_col to match the column names shown above.")
}

network_coords <- eup %>%
  group_by(Study_id, Sampling_method, Network_id) %>%
  summarise(
    lat = mean(.data[[lat_col]], na.rm = TRUE),
    lon = mean(.data[[lon_col]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(lat), !is.na(lon))

cat("\n==============================\n")
cat("Networks with coordinates:", nrow(network_coords), "\n")
cat("Studies represented:", n_distinct(network_coords$Study_id), "\n")
cat("==============================\n")

# ---- pairwise Haversine distances between networks of different studies ----
haversine_km <- function(lat1, lon1, lat2, lon2) {
  R    <- 6371
  phi1 <- lat1 * pi / 180;  phi2 <- lat2 * pi / 180
  dphi <- (lat2 - lat1) * pi / 180
  dlam <- (lon2 - lon1) * pi / 180
  a    <- sin(dphi / 2)^2 + cos(phi1) * cos(phi2) * sin(dlam / 2)^2
  2 * R * asin(sqrt(a))
}

RADIUS_KM <- 10   # adjust: smaller = stricter geographic proximity

n          <- nrow(network_coords)
pairs_list <- vector("list", n * (n - 1) / 2)
k          <- 1L

for (i in seq_len(n - 1)) {
  for (j in (i + 1):n) {
    # skip pairs from the same study or same sampling method
    if (network_coords$Study_id[i]        == network_coords$Study_id[j])        next
    if (network_coords$Sampling_method[i] == network_coords$Sampling_method[j]) next

    d <- haversine_km(network_coords$lat[i], network_coords$lon[i],
                      network_coords$lat[j], network_coords$lon[j])
    if (d <= RADIUS_KM) {
      pairs_list[[k]] <- data.frame(
        network_a  = network_coords$Network_id[i],
        study_a    = network_coords$Study_id[i],
        method_a   = network_coords$Sampling_method[i],
        network_b  = network_coords$Network_id[j],
        study_b    = network_coords$Study_id[j],
        method_b   = network_coords$Sampling_method[j],
        dist_km    = round(d, 2),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
}
close_pairs <- bind_rows(pairs_list[seq_len(k - 1)]) %>%
  arrange(dist_km)

cat(sprintf("\nNetwork pairs within %d km with different sampling methods: %d\n",
            RADIUS_KM, nrow(close_pairs)))
print(close_pairs)

# ---- method-pair summary ----
if (nrow(close_pairs) > 0) {
  method_pair_summary <- close_pairs %>%
    mutate(method_pair = paste(pmin(method_a, method_b),
                               pmax(method_a, method_b), sep = " | ")) %>%
    group_by(method_pair) %>%
    summarise(n_pairs  = n(),
              min_dist = min(dist_km),
              max_dist = max(dist_km),
              .groups  = "drop") %>%
    arrange(desc(n_pairs))

  cat("\n--- Qualifying network pairs by method combination ---\n")
  print(method_pair_summary)

  # ---- species overlap per qualifying study pair ----
  # Collapse to unique study pairs (a study may contribute multiple close networks)
  study_pairs <- close_pairs %>%
    distinct(study_a, study_b, method_a, method_b) %>%
    arrange(study_a, study_b)

  cat("\n--- Species overlap per qualifying study pair ---\n")
  for (r in seq_len(nrow(study_pairs))) {
    sa <- study_pairs$study_a[r];  sb <- study_pairs$study_b[r]
    ma <- study_pairs$method_a[r]; mb <- study_pairs$method_b[r]

    min_d <- min(close_pairs$dist_km[close_pairs$study_a == sa &
                                     close_pairs$study_b == sb])

    plants_a <- unique(na.omit(eup$Plant_accepted_name[eup$Study_id == sa]))
    plants_b <- unique(na.omit(eup$Plant_accepted_name[eup$Study_id == sb]))
    polls_a  <- unique(na.omit(eup$Pollinator_accepted_name[eup$Study_id == sa]))
    polls_b  <- unique(na.omit(eup$Pollinator_accepted_name[eup$Study_id == sb]))

    cat(sprintf("\n  %s (%s)  vs  %s (%s)  — closest networks: %.2f km\n",
                sa, ma, sb, mb, min_d))
    cat(sprintf("  Plants:      %d shared  (%d / %d)\n",
                length(intersect(plants_a, plants_b)), length(plants_a), length(plants_b)))
    cat(sprintf("  Pollinators: %d shared  (%d / %d)\n",
                length(intersect(polls_a, polls_b)), length(polls_a), length(polls_b)))
  }
}

# =============================================================================
# SECTION 2b: Detailed inspection of a specific study pair
# Set INSPECT_A and INSPECT_B to any two Study_ids from the results above.
# =============================================================================

INSPECT_A <- "40_Knight"
INSPECT_B <- "44_Knight"

inspect_pair <- function(eup, id_a, id_b) {
  ea <- eup %>% filter(Study_id == id_a)
  eb <- eup %>% filter(Study_id == id_b)
  ma <- unique(ea$Sampling_method)
  mb <- unique(eb$Sampling_method)

  cat("\n==================================================\n")
  cat(sprintf("Inspecting: %s (%s)  vs  %s (%s)\n", id_a, ma, id_b, mb))
  cat("==================================================\n")

  # ---- network inventory ----
  cat("\n--- Networks in", id_a, "---\n")
  ea %>%
    group_by(Network_id) %>%
    summarise(n_plants      = n_distinct(Plant_accepted_name),
              n_pollinators = n_distinct(Pollinator_accepted_name),
              n_links       = n(),
              lat           = mean(Latitude,  na.rm = TRUE),
              lon           = mean(Longitude, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(Network_id) %>% print()

  cat("\n--- Networks in", id_b, "---\n")
  eb %>%
    group_by(Network_id) %>%
    summarise(n_plants      = n_distinct(Plant_accepted_name),
              n_pollinators = n_distinct(Pollinator_accepted_name),
              n_links       = n(),
              lat           = mean(Latitude,  na.rm = TRUE),
              lon           = mean(Longitude, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(Network_id) %>% print()

  # ---- interaction-level overlap ----
  links_a <- ea %>%
    filter(!is.na(Plant_accepted_name), !is.na(Pollinator_accepted_name)) %>%
    distinct(plant = Plant_accepted_name, pollinator = Pollinator_accepted_name) %>%
    mutate(in_a = TRUE)

  links_b <- eb %>%
    filter(!is.na(Plant_accepted_name), !is.na(Pollinator_accepted_name)) %>%
    distinct(plant = Plant_accepted_name, pollinator = Pollinator_accepted_name) %>%
    mutate(in_b = TRUE)

  link_comparison <- full_join(links_a, links_b, by = c("plant", "pollinator")) %>%
    mutate(in_a = !is.na(in_a), in_b = !is.na(in_b),
           status = case_when(in_a & in_b ~ "shared",
                              in_a        ~ id_a,
                              TRUE        ~ id_b))

  cat("\n--- Interaction overlap (pooled across all networks) ---\n")
  link_comparison %>% count(status) %>% print()

  cat(sprintf("\n  Total unique interactions: %d\n", nrow(link_comparison)))
  cat(sprintf("  Shared by both studies:   %d  (%.1f%% of union)\n",
              sum(link_comparison$status == "shared"),
              100 * mean(link_comparison$status == "shared")))

  cat("\n--- Shared interactions ---\n")
  link_comparison %>%
    filter(status == "shared") %>%
    select(plant, pollinator) %>%
    arrange(plant, pollinator) %>%
    print(n = Inf)

  invisible(link_comparison)
}

inspect_pair(eup, INSPECT_A, INSPECT_B) # 20_Hoiss had a bunch of treatments in their study
# 11_Clough and 30_Smith - some of the networks are small, not a lot of overlap
# 40_Knight, 44_Knight - really few interactions overlap :( but the networks are large enough.
# 44_Knight is unpublished data so whether treatments were applied is unknown

# =============================================================================
# SECTION 3: FrugInt — geographic overview and candidate network identification
# Datasets: MN_1983, MN_2024, BC_seed (Pistacia only), BC_visit (Pistacia only)
# 4 candidate networks for subsequent analysis:
#   (1) MN_2024 @ Hato Ratón   (37.171, -6.332) — mist-netting
#   (2) MN_2024 @ Southern site (36.965, -6.446) — mist-netting
#   (3) BC_seed  @ Pistacia     (centroid ~1.2 km from network 2) — barcoding seeds
#   (4) BC_visit @ Pistacia     (centroid ~2.2 km from network 2) — barcoding visits
# =============================================================================

library(readr)

frugint_path <- "data/raw_data/frugint/"

# ---- load and label datasets ----
mn83 <- read_csv(paste0(frugint_path, "MN_1983.csv"), show_col_types = FALSE) %>%
  mutate(dataset = "MN_1983", method2 = "mist-netting",
         latitude = NA_real_, longitude = NA_real_,
         plantSp = plantSp, animalSp = animalSp)

mn24 <- read_csv(paste0(frugint_path, "MN_2024.csv"), show_col_types = FALSE) %>%
  filter(grepl("Pistacia", vegetation)) %>%
  mutate(dataset = "MN_2024",
         network = ifelse(latitude > 37.1, "MN_2024_HatoRaton", "MN_2024_South"))

bcs <- read_csv(paste0(frugint_path, "BC_seed.csv"), show_col_types = FALSE) %>%
  filter(grepl("Pistacia", vegetation)) %>%
  mutate(dataset  = "BC_seed",
         network  = "BC_seed_Pistacia",
         latitude = as.numeric(latitude),
         longitude = as.numeric(longitude))

bcv <- read_csv(paste0(frugint_path, "BC_visit.csv"), show_col_types = FALSE) %>%
  filter(grepl("Pistacia", vegetation)) %>%
  mutate(dataset  = "BC_visit",
         network  = "BC_visit_Pistacia",
         latitude = as.numeric(latitude),
         longitude = as.numeric(longitude))

# ---- text: site summary ----
haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371; d2r <- pi / 180
  dlat <- (lat2 - lat1) * d2r; dlon <- (lon2 - lon1) * d2r
  a <- sin(dlat/2)^2 + cos(lat1*d2r)*cos(lat2*d2r)*sin(dlon/2)^2
  2 * R * asin(sqrt(a))
}

sites <- bind_rows(
  mn24 %>% distinct(dataset, network, method2, latitude, longitude),
  bcs  %>% group_by(dataset, network, method2) %>%
           summarise(latitude = mean(latitude, na.rm = TRUE),
                     longitude = mean(longitude, na.rm = TRUE), .groups = "drop"),
  bcv  %>% group_by(dataset, network, method2) %>%
           summarise(latitude = mean(latitude, na.rm = TRUE),
                     longitude = mean(longitude, na.rm = TRUE), .groups = "drop")
)

cat("\n==============================\n")
cat("SECTION 3: FrugInt candidate networks\n")
cat("==============================\n")
cat("\n--- Site centroids ---\n")
sites %>% mutate(across(c(latitude, longitude), ~round(.x, 5))) %>% print()

cat("\n--- Pairwise distances between candidate networks (km) ---\n")
for (i in seq_len(nrow(sites) - 1)) {
  for (j in (i + 1):nrow(sites)) {
    d <- haversine_km(sites$latitude[i], sites$longitude[i],
                      sites$latitude[j], sites$longitude[j])
    cat(sprintf("  %s  <->  %s : %.2f km\n",
                sites$network[i], sites$network[j], round(d, 2)))
  }
}

# ---- text: network size per candidate ----
cat("\n--- Network size per candidate (plants × animals × interactions) ---\n")

network_size <- function(df, net_col = "network") {
  df %>%
    group_by(network = .data[[net_col]]) %>%
    summarise(n_plants  = n_distinct(plantSp),
              n_animals = n_distinct(animalSp),
              n_links   = n_distinct(paste(plantSp, animalSp)),
              .groups   = "drop")
}

bind_rows(
  network_size(mn24),
  network_size(bcs),
  network_size(bcv)
) %>% print()

# ---- text: species overlap between candidate networks ----
get_sp <- function(df, col) unique(na.omit(df[[col]]))

nets <- list(
  MN_2024_HatoRaton  = mn24 %>% filter(network == "MN_2024_HatoRaton"),
  MN_2024_South      = mn24 %>% filter(network == "MN_2024_South"),
  BC_seed_Pistacia   = bcs,
  BC_visit_Pistacia  = bcv
)

cat("\n--- Plant species overlap between candidate networks ---\n")
nms <- names(nets)
for (i in seq_along(nms)) {
  for (j in seq_along(nms)) {
    if (j <= i) next
    p1 <- get_sp(nets[[i]], "plantSp");  p2 <- get_sp(nets[[j]], "plantSp")
    cat(sprintf("  %s ∩ %s : %d shared plants (%d / %d)\n",
                nms[i], nms[j], length(intersect(p1, p2)), length(p1), length(p2)))
  }
}

cat("\n--- Animal species overlap between candidate networks ---\n")
for (i in seq_along(nms)) {
  for (j in seq_along(nms)) {
    if (j <= i) next
    a1 <- get_sp(nets[[i]], "animalSp"); a2 <- get_sp(nets[[j]], "animalSp")
    cat(sprintf("  %s ∩ %s : %d shared animals (%d / %d)\n",
                nms[i], nms[j], length(intersect(a1, a2)), length(a1), length(a2)))
  }
}

# ---- map: all sampling points coloured by dataset, shaped by method2 ----
map_df <- bind_rows(
  mn24 %>% select(dataset, method2, network, latitude, longitude),
  bcs  %>% select(dataset, method2, network, latitude, longitude),
  bcv  %>% select(dataset, method2, network, latitude, longitude)
) %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  mutate(
    dataset  = factor(dataset,  levels = c("MN_2024", "BC_seed", "BC_visit")),
    method2  = factor(method2,  levels = c("mist-netting", "seed-underPlant", "visit"))
  )

# site-level centroids for labels
site_labels <- map_df %>%
  group_by(network, dataset, method2) %>%
  summarise(lat = mean(latitude), lon = mean(longitude), .groups = "drop")

spain <- map_data("world", region = "Spain")

map_plot <- ggplot() +
  geom_polygon(data = spain, aes(x = long, y = lat, group = group),
               fill = "grey90", colour = "grey60", linewidth = 0.3) +
  geom_point(data = map_df,
             aes(x = longitude, y = latitude, colour = dataset, shape = method2),
             alpha = 0.5, size = 1.5) +
  geom_point(data = site_labels,
             aes(x = lon, y = lat, colour = dataset, shape = method2),
             size = 4, stroke = 1.2) +
  ggrepel::geom_label_repel(data = site_labels,
             aes(x = lon, y = lat, label = network, colour = dataset),
             size = 3, show.legend = FALSE, box.padding = 0.4) +
  coord_fixed(ratio = 1,
              xlim = c(min(map_df$longitude) - 0.05, max(map_df$longitude) + 0.05),
              ylim = c(min(map_df$latitude)  - 0.05, max(map_df$latitude)  + 0.10)) +
  scale_colour_manual(values = c("MN_2024"  = "#e63946",
                                 "BC_seed"  = "#2a9d8f",
                                 "BC_visit" = "#f4a261")) +
  scale_shape_manual(values = c("mist-netting"    = 16,
                                "seed-underPlant" = 17,
                                "visit"           = 15)) +
  labs(title   = "FrugInt candidate networks — Pistacia-dominated scrubland",
       x = "Longitude", y = "Latitude",
       colour = "Dataset", shape = "Method") +
  theme_minimal(base_size = 11) +
  theme(legend.position  = "right",
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

# use ggrepel if available, otherwise fall back to geom_text
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  map_plot <- ggplot() +
    geom_polygon(data = spain, aes(x = long, y = lat, group = group),
                 fill = "grey90", colour = "grey60", linewidth = 0.3) +
    geom_point(data = map_df,
               aes(x = longitude, y = latitude, colour = dataset, shape = method2),
               alpha = 0.5, size = 1.5) +
    geom_point(data = site_labels,
               aes(x = lon, y = lat, colour = dataset, shape = method2),
               size = 4, stroke = 1.2) +
    geom_text(data = site_labels,
              aes(x = lon, y = lat, label = network, colour = dataset),
              vjust = -1.2, size = 3, show.legend = FALSE) +
    coord_fixed(ratio = 1,
                xlim = c(min(map_df$longitude) - 0.05, max(map_df$longitude) + 0.05),
                ylim = c(min(map_df$latitude)  - 0.05, max(map_df$latitude)  + 0.10)) +
    scale_colour_manual(values = c("MN_2024"  = "#e63946",
                                   "BC_seed"  = "#2a9d8f",
                                   "BC_visit" = "#f4a261")) +
    scale_shape_manual(values = c("mist-netting"    = 16,
                                  "seed-underPlant" = 17,
                                  "visit"           = 15)) +
    labs(title   = "FrugInt candidate networks — Pistacia-dominated scrubland",
         x = "Longitude", y = "Latitude",
         colour = "Dataset", shape = "Method") +
    theme_minimal(base_size = 11) +
    theme(legend.position  = "right",
          panel.grid.minor = element_blank(),
          plot.title       = element_text(face = "bold"))
}

print(map_plot)
