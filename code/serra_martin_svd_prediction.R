# ---- Serra-Martin SVD link prediction ----
# Predict plant-pollinator interactions for each of 6 habitat patches
# using SVD-based matrix completion (softImpute). For each focal site,
# the aggregate of all other patches under the same sampling method is used
# for training. Outputs one prediction CSV per site × method for downstream
# classification in serra_martin_link_classification.R.

source("code/common.R")
library(tidyverse)
library(softImpute)
library(pROC)
library(PRROC)

## ---- 1. Parameters ----
n_sim <- 50   # bootstrap iterations
set.seed(42)

## ---- 2. Load and prepare Serra-Martin data ----
df_raw <- read_delim(
  "data/raw_data/serra_martin_pollination/cabrera_22_23_habitat.csv",
  delim = ";", show_col_types = FALSE
)

df_clean <- df_raw %>%
  filter(!is.na(visita)) %>%                          # drop 50 blank rows
  mutate(`Plant sp` = recode(`Plant sp`,
    "Daucus carota L. subsp. Majoricus" = "Daucus carota",
    "Rosmarinus officinalis"            = "Salvia rosmarinus"
  )) %>%                                              # harmonise plant names
  filter(!`Plant sp` %in% c("Anacamptis pyramidalis", "Gladiolus communis")) %>%  # remove camera-training orchids
  filter(!is.na(Pollinator), Pollinator != "")        # keep only rows with a recorded pollinator

sites   <- sort(unique(df_clean$habitat))
methods <- sort(unique(df_clean$Method))

# Format as emln-style long-format data (one data frame per method).
# layer_from = habitat patch, node_from = plant, node_to = pollinator, weight = 1.
# Deduplicated to unique (site, plant, pollinator) combinations — binary presence.
make_aggregated_df <- function(data, method) {
  data %>%
    filter(Method == method) %>%
    distinct(habitat, `Plant sp`, Pollinator) %>%
    transmute(
      layer_from = habitat,
      node_from  = `Plant sp`,
      layer_to   = habitat,
      node_to    = Pollinator,
      weight     = 1L
    )
}

## ---- 3. Prediction ----
results_file <- "results/serra_martin_predictions.rds"

combined_results <- data.frame()

if (file.exists(results_file)) {
  print("Existing results file found — reading the file and proceeding to analysis")
  combined_results <- readRDS(results_file)
  print("finished loading prediction results")

} else {
  for (m in methods) {
    print(paste("=== Method:", m, "==="))
    aggregated_df <- make_aggregated_df(df_clean, m)

    for (focal_site in sites) {
      other_sites <- setdiff(sites, focal_site)
      print(paste("** training:", paste(other_sites, collapse = ", "),
                  " | predicting:", focal_site, "**"))

      # Build training matrix A (aggregate of all other sites under this method)
      A <- build_interaction_matrix(data = aggregated_df, layers_to_filter = other_sites)

      # Build prediction matrix P (focal site only)
      P <- build_interaction_matrix(data = aggregated_df, layers_to_filter = focal_site)

      node_to   <- rownames(P)
      node_from <- colnames(P)

      # Map 0s and 1s in P
      # max(1, ...) ensures at least one link is withheld even for sparse rpi sites
      num_1_to_remove <- max(1L, floor(sum(P > 0, na.rm = TRUE) * prop_ones_to_remove))
      ones_in_P  <- which(P > 0, arr.ind = TRUE)

      num_0_to_remove <- num_1_to_remove
      prop_0_removed  <- num_0_to_remove / sum(P == 0, na.rm = TRUE)
      zeros_in_P <- which(P == 0, arr.ind = TRUE)

      print(paste("1 remove:", num_1_to_remove))
      print(paste("all 1   :", nrow(ones_in_P)))
      print(paste("0s to remove:", num_0_to_remove))
      print(paste("all zeros   :", nrow(zeros_in_P)))
      print(paste("prop of zeros removed:", round(prop_0_removed, 3)))

      bootstrapping_results <- NULL
      P_original <- P

      for (i in 1:n_sim) {
        # Remove 1s
        remove_indices          <- ones_in_P[sample(1:nrow(ones_in_P),  num_1_to_remove), ]
        P[remove_indices]       <- NA

        # Sample 0s
        zeros_to_remove_indices   <- zeros_in_P[sample(1:nrow(zeros_in_P), num_0_to_remove), ]
        P[zeros_to_remove_indices] <- NA

        ### ---- Build combined matrix C ----
        all_row_ids <- unique(c(rownames(A), rownames(P)))
        all_col_ids <- unique(c(colnames(A), colnames(P)))
        C <- matrix(0, nrow = length(all_row_ids), ncol = length(all_col_ids),
                    dimnames = list(all_row_ids, all_col_ids))

        C[rownames(A), colnames(A)] <- A

        # Training and prediction always differ (leave-one-out): sum overlapping entries
        C[rownames(P), colnames(P)] <- ifelse(
          is.na(C[rownames(P), colnames(P)]),
          NA,
          C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)]
        )

        # biScale centering
        C <- biScale(C, row.center = TRUE, col.center = TRUE,
                     row.scale = FALSE, col.scale = FALSE)
        row_centers <- attr(C, "biScale:row")$center
        col_centers <- attr(C, "biScale:column")$center

        ### ---- Predict with SVD for all k/lambda combinations ----
        k_values      <- c(2, 5, 10)
        lam0          <- lambda0(C)
        lambda_values <- c(1, 5, 50, 100, lam0)

        results <- data.frame(k = integer(), lambda = numeric(),
                              original_links = numeric(), predicted_values = numeric(),
                              input_lambda = numeric())
        not_removed_all <- NULL

        for (k in k_values) {
          for (lambda in lambda_values) {
            r <- implement_impute(C, k, lambda,
                                  P, remove_indices, zeros_to_remove_indices, P_original,
                                  back_trans_values = list(row_centers = row_centers,
                                                           col_centers = col_centers))
            r$results$input_lambda     <- lambda
            r$not_removed$input_lambda <- lambda
            results         <- rbind(results, r$results)
            not_removed_all <- rbind(not_removed_all, r$not_removed)
          }
        }

        complete_edges_all     <- rbind(results, not_removed_all)
        complete_edges_all$itr <- i
        bootstrapping_results  <- rbind(bootstrapping_results, complete_edges_all)

        P <- P_original
      }

      combined_results <- rbind(
        combined_results,
        cbind(
          data.frame(
            method              = m,
            focal_site          = focal_site,
            prop_ones_removed   = prop_ones_to_remove,
            amount_of_removed_1 = num_1_to_remove,
            amount_of_removed_0 = num_0_to_remove,
            prop_0_removed      = prop_0_removed
          ),
          bootstrapping_results
        )
      )
    }
  }
  saveRDS(combined_results, file = results_file)
}

# Filter to optimal k/lambda (k=2, lambda=lambda0 only)
combined_results <- combined_results %>%
  filter(k == 2) %>%
  filter(!(input_lambda %in% c(1, 5, 50, 100)))

## ---- 4. Analysis ----
df <- combined_results %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

### ---- Select optimal threshold (max F0.5) ----
thresholds <- seq(0, 1, by = 0.1)

df_prepped <- df %>%
  filter(removed == 1) %>%
  mutate(
    predicted_prob  = sigmoid(predicted_values),
    original_binary = if_else(original_links > 0, 1, 0)
  )

df_thresh <- df_prepped %>%
  tidyr::expand_grid(threshold = thresholds) %>%
  mutate(predicted_bin = if_else(predicted_prob > threshold, 1, 0)) %>%
  group_by(method, focal_site, itr, threshold) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin == 1),
    FN = sum(original_binary == 1 & predicted_bin == 0),
    TN = sum(original_binary == 0 & predicted_bin == 0),
    FP = sum(original_binary == 0 & predicted_bin == 1),
    specificity       = TN / (TN + FP),
    precision         = TP / (TP + FP),
    recall            = TP / (TP + FN),
    f05_score         = (1.25) * (precision * recall) / ((0.25 * precision) + recall),
    balanced_accuracy = (recall + specificity) / 2,
    mcc  = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
    mse  = mean((predicted_values - original_links)^2, na.rm = TRUE),
    rmse = sqrt(mse)
  ) %>%
  ungroup() %>%
  group_by(method, focal_site, threshold) %>%
  summarise(
    TP = mean(TP, na.rm = TRUE), FN = mean(FN, na.rm = TRUE),
    TN = mean(TN, na.rm = TRUE), FP = mean(FP, na.rm = TRUE),
    specificity       = mean(specificity,       na.rm = TRUE),
    precision         = mean(precision,         na.rm = TRUE),
    recall            = mean(recall,            na.rm = TRUE),
    f05_score         = mean(f05_score,         na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    mcc               = mean(mcc,               na.rm = TRUE),
    mse               = mean(mse,               na.rm = TRUE),
    rmse              = mean(rmse,              na.rm = TRUE)
  ) %>%
  ungroup()

df_avg <- df_thresh %>%
  group_by(threshold) %>%
  summarise(across(
    c(specificity, precision, recall, f05_score, balanced_accuracy, mcc),
    mean, na.rm = TRUE
  )) %>%
  pivot_longer(-threshold, names_to = "metric", values_to = "value")

df_wide <- df_avg %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  arrange(threshold)

best_discrete <- df_wide %>% slice_max(f05_score, n = 1)
best_discrete_threshold <- best_discrete$threshold
best_discrete_threshold

df_removed <- df %>%
  filter(removed == 1) %>%
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%
  mutate(predicted_bin_sigm  = if_else(predicted_prob_sigm > best_discrete_threshold, 1, 0)) %>%
  mutate(original_binary     = if_else(original_links > 0, 1, 0))

result_summary <- df_removed %>%
  group_by(method, focal_site, itr) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin_sigm == 1),
    FN = sum(original_binary == 1 & predicted_bin_sigm == 0),
    TN = sum(original_binary == 0 & predicted_bin_sigm == 0),
    FP = sum(original_binary == 0 & predicted_bin_sigm == 1),
    specificity       = TN / (TN + FP),
    precision         = TP / (TP + FP),
    recall            = TP / (TP + FN),
    f05_score         = (1.25) * (precision * recall) / ((0.25 * precision) + recall),
    balanced_accuracy = (recall + specificity) / 2,
    nse  = 1 - sum((predicted_values - original_links)^2, na.rm = TRUE) /
      sum((original_links - mean(original_links, na.rm = TRUE))^2, na.rm = TRUE),
    nnse = 1 / (2 - nse)
  ) %>%
  ungroup() %>%
  group_by(method, focal_site) %>%
  summarise(
    TP = mean(TP, na.rm = TRUE), FN = mean(FN, na.rm = TRUE),
    TN = mean(TN, na.rm = TRUE), FP = mean(FP, na.rm = TRUE),
    specificity       = mean(specificity,       na.rm = TRUE),
    precision         = mean(precision,         na.rm = TRUE),
    recall            = mean(recall,            na.rm = TRUE),
    f05_score         = mean(f05_score,         na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    nse               = mean(nse,               na.rm = TRUE),
    nnse              = mean(nnse,              na.rm = TRUE)
  ) %>%
  ungroup()

head(result_summary)
summary(result_summary)

## ---- 5. Final full prediction (no link withholding) ----
# A single SVD pass per site × method to get predictions for ALL pairs in P,
# including zeros never sampled during bootstrapping. Used for saving CSVs only.

final_predictions_file <- "results/serra_martin_final_predictions.rds"

all_final_preds <- data.frame()

if (file.exists(final_predictions_file)) {
  print("Existing final predictions found — loading from disk")
  all_final_preds <- readRDS(final_predictions_file)
} else {
  for (m in methods) {
    aggregated_df <- make_aggregated_df(df_clean, m)

    for (focal_site in sites) {
      other_sites <- setdiff(sites, focal_site)

      A      <- build_interaction_matrix(data = aggregated_df, layers_to_filter = other_sites)
      P_full <- build_interaction_matrix(data = aggregated_df, layers_to_filter = focal_site)

      all_row_ids <- unique(c(rownames(A), rownames(P_full)))
      all_col_ids <- unique(c(colnames(A), colnames(P_full)))
      C <- matrix(0, nrow = length(all_row_ids), ncol = length(all_col_ids),
                  dimnames = list(all_row_ids, all_col_ids))
      C[rownames(A),      colnames(A)]      <- A
      C[rownames(P_full), colnames(P_full)] <- C[rownames(P_full), colnames(P_full)] + P_full

      C <- biScale(C, row.center = TRUE, col.center = TRUE,
                   row.scale = FALSE, col.scale = FALSE)
      row_centers <- attr(C, "biScale:row")$center
      col_centers <- attr(C, "biScale:column")$center

      lam0 <- lambda0(C)
      fit  <- softImpute(C, rank.max = 2, lambda = lam0, type = "svd", maxit = 600)
      C_rec      <- softImpute::complete(C, fit)
      C_rec_orig <- C_rec + outer(row_centers[rownames(C)], col_centers[colnames(C)], "+")
      P_rec      <- C_rec_orig[rownames(P_full), colnames(P_full)]

      site_preds <- expand.grid(
        node_to   = rownames(P_full),
        node_from = colnames(P_full),
        stringsAsFactors = FALSE
      ) %>%
        mutate(
          method           = m,
          focal_site       = focal_site,
          original_links   = P_full[cbind(node_to, node_from)],
          predicted_values = P_rec[cbind(node_to, node_from)]
        )

      all_final_preds <- rbind(all_final_preds, site_preds)
    }
  }
  saveRDS(all_final_preds, final_predictions_file)
}

## ---- 6. Save per-site × method prediction CSVs ----
dir.create("results/predictions", showWarnings = FALSE, recursive = TRUE)

all_final_preds <- all_final_preds %>%
  mutate(
    predicted_values    = if_else(predicted_values < 0, 0, predicted_values),
    predicted_prob_sigm = sigmoid(predicted_values),
    predicted_bin_sigm  = if_else(predicted_prob_sigm > best_discrete_threshold, 1L, 0L),
    original_binary     = if_else(original_links > 0, 1L, 0L)
  )

for (m in methods) {
  for (site in sites) {
    df_site <- all_final_preds %>%
      filter(method == m, focal_site == site)

    if (nrow(df_site) == 0) next

    df_out <- df_site %>%
      transmute(
        link_ID        = row_number() - 1L,
        higher_level   = node_to,            # pollinator
        lower_level    = node_from,          # plant
        ground_truth   = original_binary,
        probability    = predicted_prob_sigm,
        classification = predicted_bin_sigm
      )

    fname <- sprintf(
      "results/predictions/serra_martin_%s_%s_prediction_results.csv",
      gsub(" ", "_", tolower(site)), m
    )
    write_csv(df_out, fname)
    print(paste("Saved:", fname))
  }
}
