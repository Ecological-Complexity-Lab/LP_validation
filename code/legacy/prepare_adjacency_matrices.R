## ---- Prepare adjacency matrices for link prediction ----
# Reads MB-LDB and MB-RDB sheets from DryadMetabarcodingData.xlsx and produces
# one binary adjacency matrix (plants x pollinators) per location x method
# combination, plus combined metaweb matrices.
# Methods: "observation" (visual, Observation column) — same in both sheets,
#           saved once without a database suffix.
#          "metabarcoding" (pollen DNA, Taxon1-Taxon9) — saved separately for
#           LDB and RDB with suffixes _LDB / _RDB.
# Months and sexes are collapsed within each location (distinct pairs only).
# Output: CSV files in data/adjacency_matrices/

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

### ---- pollination ----
### ---- paths ----
raw_path <- "data/raw_data/DryadMetabarcodingData.xlsx"
out_dir  <- "data/adjacency_matrices"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

### ---- get distinct plant-pollinator pairs from a long-format data frame ----
get_pairs <- function(df, plant_col, poll_col) {
  df %>%
    select(plant = all_of(plant_col), pollinator = all_of(poll_col)) %>%
    filter(!is.na(plant), !is.na(pollinator), plant != "", pollinator != "") %>%
    distinct()
}

### ---- build binary adjacency matrix from distinct pairs ----
build_adj <- function(pairs) {
  plants      <- sort(unique(pairs$plant))
  pollinators <- sort(unique(pairs$pollinator))

  mat <- matrix(0L,
                nrow = length(plants),
                ncol = length(pollinators),
                dimnames = list(plants, pollinators))

  for (i in seq_len(nrow(pairs))) {
    mat[pairs$plant[i], pairs$pollinator[i]] <- 1L
  }

  as.data.frame(mat)
}

### ---- sanity checks ----
sanity_check <- function(mat, pairs, csv_path, label) {
  pass <- TRUE
  pfx  <- paste0("[", label, "]")

  check <- function(ok, msg) {
    cat(pfx, if (ok) "PASS" else "FAIL", "--", msg, "\n")
    if (!ok) pass <<- FALSE
  }

  # 1. link count: number of 1s must equal number of distinct pairs
  check(sum(mat) == nrow(pairs),
        paste("link count: matrix has", sum(mat), "ones,",
              nrow(pairs), "distinct pairs in raw data"))

  # 2. all values are 0 or 1
  check(all(as.matrix(mat) %in% c(0L, 1L)),
        "all matrix values are 0 or 1")

  # 3. row names cover exactly the plants in the raw pairs
  check(setequal(rownames(mat), unique(pairs$plant)),
        paste("plant species coverage:", nrow(mat), "rows"))

  # 4. column names cover exactly the pollinators in the raw pairs
  check(setequal(colnames(mat), unique(pairs$pollinator)),
        paste("pollinator species coverage:", ncol(mat), "columns"))

  # 5. no all-zero row (every plant has at least one interaction)
  empty_rows <- rownames(mat)[rowSums(mat) == 0]
  check(length(empty_rows) == 0,
        paste("no empty rows (plants with zero links);",
              if (length(empty_rows)) paste("empty:", paste(empty_rows, collapse = ", ")) else ""))

  # 6. no all-zero column (every pollinator has at least one interaction)
  empty_cols <- colnames(mat)[colSums(mat) == 0]
  check(length(empty_cols) == 0,
        paste("no empty columns (pollinators with zero links);",
              if (length(empty_cols)) paste("empty:", paste(empty_cols, collapse = ", ")) else ""))

  # 7. spot-check: 5 random raw pairs must be 1 in the matrix; print each example
  set.seed(42)
  sample_idx <- sample(seq_len(nrow(pairs)), min(5, nrow(pairs)))
  spot_results <- vapply(sample_idx, function(i) {
    val <- mat[pairs$plant[i], pairs$pollinator[i]]
    cat(pfx, "  example:", pairs$pollinator[i], "<->", pairs$plant[i],
        "=> matrix value:", val, "\n")
    val == 1L
  }, logical(1))
  check(all(spot_results), "spot-check: all 5 sampled raw pairs are 1 in matrix")

  # 8. round-trip: re-read CSV, dimensions and link count must match
  rt <- read.csv(csv_path, row.names = 1, check.names = FALSE)
  check(nrow(rt) == nrow(mat) && ncol(rt) == ncol(mat),
        paste("round-trip CSV dimensions match:", nrow(rt), "x", ncol(rt)))
  check(sum(rt) == sum(mat),
        paste("round-trip CSV link count matches:", sum(rt)))

  if (!pass) stop("Sanity checks failed for ", label, ". See messages above.")
  invisible(NULL)
}

### ---- load data ----
dat_ldb <- read_excel(raw_path, sheet = "MB-LDB") %>%
  mutate(pollinator = str_trim(paste(Genus, Species)))

dat_rdb <- read_excel(raw_path, sheet = "MB-RDB") %>%
  mutate(pollinator = str_trim(paste(Genus, Species)))

locations <- unique(dat_ldb$Location)

### ---- method 1: visual observation (identical in both sheets — save once) ----
obs_long <- dat_ldb %>%
  select(pollinator, Location, plant = Observation)

for (loc in locations) {
  pairs_obs <- get_pairs(filter(obs_long, Location == loc), "plant", "pollinator")
  mat_obs   <- build_adj(pairs_obs)
  csv_obs   <- file.path(out_dir, paste0(loc, "_observation.csv"))
  write.csv(mat_obs, file = csv_obs, row.names = TRUE)
  sanity_check(mat_obs, pairs_obs, csv_obs, paste(loc, "observation"))
  cat("Saved", loc, "observation |",
      nrow(mat_obs), "plants x", ncol(mat_obs), "pollinators\n\n")
}

# metaweb observation
pairs_meta_obs <- get_pairs(obs_long, "plant", "pollinator")
mat_meta_obs   <- build_adj(pairs_meta_obs)
csv_meta_obs   <- file.path(out_dir, "metaweb_observation.csv")
write.csv(mat_meta_obs, file = csv_meta_obs, row.names = TRUE)
sanity_check(mat_meta_obs, pairs_meta_obs, csv_meta_obs, "metaweb observation")
cat("Saved metaweb observation |",
    nrow(mat_meta_obs), "plants x", ncol(mat_meta_obs), "pollinators\n\n")

### ---- method 2: metabarcoding (LDB and RDB) ----
databases <- list(LDB = dat_ldb, RDB = dat_rdb)

for (db_name in names(databases)) {
  dat_db <- databases[[db_name]]

  mb_long <- dat_db %>%
    select(pollinator, Location,
           Taxon1, Taxon2, Taxon3, Taxon4, Taxon5,
           Taxon6, Taxon7, Taxon8, Taxon9) %>%
    pivot_longer(cols      = starts_with("Taxon"),
                 names_to  = "taxon_col",
                 values_to = "plant") %>%
    select(pollinator, Location, plant)

  for (loc in locations) {
    pairs_mb <- get_pairs(filter(mb_long, Location == loc), "plant", "pollinator")
    mat_mb   <- build_adj(pairs_mb)
    csv_mb   <- file.path(out_dir, paste0(loc, "_metabarcoding_", db_name, ".csv"))
    write.csv(mat_mb, file = csv_mb, row.names = TRUE)
    sanity_check(mat_mb, pairs_mb, csv_mb, paste(loc, "metabarcoding", db_name))
    cat("Saved", loc, "metabarcoding", db_name, "|",
        nrow(mat_mb), "plants x", ncol(mat_mb), "pollinators\n\n")
  }

  # metaweb metabarcoding
  pairs_meta_mb <- get_pairs(mb_long, "plant", "pollinator")
  mat_meta_mb   <- build_adj(pairs_meta_mb)
  csv_meta_mb   <- file.path(out_dir, paste0("metaweb_metabarcoding_", db_name, ".csv"))
  write.csv(mat_meta_mb, file = csv_meta_mb, row.names = TRUE)
  sanity_check(mat_meta_mb, pairs_meta_mb, csv_meta_mb,
               paste("metaweb metabarcoding", db_name))
  cat("Saved metaweb metabarcoding", db_name, "|",
      nrow(mat_meta_mb), "plants x", ncol(mat_meta_mb), "pollinators\n\n")
}

cat("Done. Files written to:", out_dir, "\n")

# ---- Frugint ----
# FrugInt frugivory networks — binary adjacency matrices (plants × animals)
# Networks:
#   frugint_MN2024_HatoRaton  — mist-netting, lat > 37.1 (Hato Ratón)
#   frugint_MN2024_South      — mist-netting, lat < 37.1 (southern Pistacia site)
#   frugint_MN2024_metaweb    — both MN_2024 Pistacia sites combined
#   frugint_BCvisit_Pistacia  — barcoding visits, Pistacia habitat (all pooled)
#   frugint_BCseed_Pistacia   — barcoding seeds, Pistacia habitat (all method2 pooled)

library(readr)

frugint_path <- "data/raw_data/frugint/"

mn24 <- read_csv(paste0(frugint_path, "MN_2024.csv"), show_col_types = FALSE) %>%
  filter(grepl("Pistacia", vegetation)) %>%
  mutate(latitude = as.numeric(latitude))

bcv <- read_csv(paste0(frugint_path, "BC_visit.csv"), show_col_types = FALSE) %>%
  filter(grepl("Pistacia", vegetation))

bcs <- read_csv(paste0(frugint_path, "BC_seed.csv"), show_col_types = FALSE) %>%
  filter(grepl("Pistacia", vegetation))

# frugint uses plantSp / animalSp — reuse existing get_pairs / build_adj / sanity_check
frugint_networks <- list(
  frugint_MN2024_HatoRaton = mn24 %>% filter(latitude > 37.1),
  frugint_MN2024_South     = mn24 %>% filter(latitude < 37.1),
  frugint_MN2024_metaweb   = mn24,
  frugint_BCvisit_Pistacia = bcv,
  frugint_BCseed_Pistacia  = bcs
)

for (net_name in names(frugint_networks)) {
  df_net <- frugint_networks[[net_name]]

  pairs <- get_pairs(df_net, "plantSp", "animalSp")

  if (nrow(pairs) == 0) {
    cat("SKIP", net_name, "— no interactions after filtering\n\n")
    next
  }

  mat     <- build_adj(pairs)
  csv_out <- file.path(out_dir, paste0(net_name, ".csv"))
  write.csv(mat, file = csv_out, row.names = TRUE)
  sanity_check(mat, pairs, csv_out, net_name)
  cat("Saved", net_name, "|",
      nrow(mat), "plants x", ncol(mat), "animals\n\n")
}

cat("FrugInt matrices written to:", out_dir, "\n")
