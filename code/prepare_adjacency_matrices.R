## ---- Prepare adjacency matrices for link prediction ----
# Reads MB-LDB sheet from DryadMetabarcodingData.xlsx and produces one binary
# adjacency matrix (plants x pollinators) per location x method combination.
# Methods: "observation" (visual, Observation column) and
#          "metabarcoding" (pollen DNA, Taxon1-Taxon9 columns, LDB database).
# Months and sexes are collapsed within each location (distinct pairs only).
# Output: CSV files in results/adjacency_matrices/

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

# ---- paths ----
raw_path <- "data/raw_data/DryadMetabarcodingData.xlsx"
out_dir  <- "data/adjacency_matrices"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- load data ----
dat <- read_excel(raw_path, sheet = "MB-LDB")

# pollinator label: Genus + Species only (sex and month collapsed downstream)
dat <- dat %>%
  mutate(pollinator = str_trim(paste(Genus, Species)))

locations <- unique(dat$Location)

# ---- get distinct plant-pollinator pairs from a long-format data frame ----
get_pairs <- function(df, plant_col, poll_col) {
  df %>%
    select(plant = all_of(plant_col), pollinator = all_of(poll_col)) %>%
    filter(!is.na(plant), !is.na(pollinator), plant != "", pollinator != "") %>%
    distinct()
}

# ---- build binary adjacency matrix from distinct pairs ----
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

# ---- sanity checks ----
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

# ---- method 1: visual observation ----
obs_long <- dat %>%
  select(pollinator, Location, plant = Observation)

# ---- method 2: metabarcoding (LDB) ----
mb_long <- dat %>%
  select(pollinator, Location,
         Taxon1, Taxon2, Taxon3, Taxon4, Taxon5,
         Taxon6, Taxon7, Taxon8, Taxon9) %>%
  pivot_longer(cols      = starts_with("Taxon"),
               names_to  = "taxon_col",
               values_to = "plant") %>%
  select(pollinator, Location, plant)

# ---- build, save, and check one CSV per location x method ----
for (loc in locations) {

  pairs_obs <- get_pairs(filter(obs_long, Location == loc), "plant", "pollinator")
  mat_obs   <- build_adj(pairs_obs)
  csv_obs   <- file.path(out_dir, paste0(loc, "_observation.csv"))
  write.csv(mat_obs, file = csv_obs, row.names = TRUE)
  sanity_check(mat_obs, pairs_obs, csv_obs, paste(loc, "observation"))

  pairs_mb <- get_pairs(filter(mb_long, Location == loc), "plant", "pollinator")
  mat_mb   <- build_adj(pairs_mb)
  csv_mb   <- file.path(out_dir, paste0(loc, "_metabarcoding.csv"))
  write.csv(mat_mb, file = csv_mb, row.names = TRUE)
  sanity_check(mat_mb, pairs_mb, csv_mb, paste(loc, "metabarcoding"))

  cat("Saved", loc,
      "| obs:", nrow(mat_obs), "plants x", ncol(mat_obs), "pollinators",
      "| mb:",  nrow(mat_mb),  "plants x", ncol(mat_mb),  "pollinators\n\n")
}

# ---- metaweb: all locations combined, one matrix per method ----
pairs_meta_obs <- get_pairs(obs_long, "plant", "pollinator")
mat_meta_obs   <- build_adj(pairs_meta_obs)
csv_meta_obs   <- file.path(out_dir, "metaweb_observation.csv")
write.csv(mat_meta_obs, file = csv_meta_obs, row.names = TRUE)
sanity_check(mat_meta_obs, pairs_meta_obs, csv_meta_obs, "metaweb observation")

pairs_meta_mb <- get_pairs(mb_long, "plant", "pollinator")
mat_meta_mb   <- build_adj(pairs_meta_mb)
csv_meta_mb   <- file.path(out_dir, "metaweb_metabarcoding.csv")
write.csv(mat_meta_mb, file = csv_meta_mb, row.names = TRUE)
sanity_check(mat_meta_mb, pairs_meta_mb, csv_meta_mb, "metaweb metabarcoding")

cat("Saved metaweb",
    "| obs:", nrow(mat_meta_obs), "plants x", ncol(mat_meta_obs), "pollinators",
    "| mb:",  nrow(mat_meta_mb),  "plants x", ncol(mat_meta_mb),  "pollinators\n\n")

cat("Done. Files written to:", out_dir, "\n")
