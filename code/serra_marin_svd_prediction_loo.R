# ---- Serra-Marin SVD self-prediction (leave-one-out) ----
# For each site × method, predict every interaction in P using only the
# within-site matrix structure — no cross-site training. Every entry in P is
# withheld once, predicted by SVD from the remaining entries, then restored.
# Equivalent to the self-prediction (/2) branch of svd_link_prediction.R.
# Outputs a single combined CSV across all sites and methods.

source("code/common.R")
library(tidyverse)
library(softImpute)
library(pROC)
library(PRROC)

## ---- 1. Parameters ----
set.seed(42)

## ---- 2. Load and prepare Serra-Marin data ----
df_raw <- read_delim(
  "data/raw_data/serra_marin_pollination/cabrera_22_23_habitat.csv",
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

## ---- 3. Leave-one-out self-prediction ----
results_file <- "results/serra_marin_loo_predictions.rds"

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
      print(paste("  Site:", focal_site))

      P          <- build_interaction_matrix(data = aggregated_df, layers_to_filter = focal_site)
      P_original <- P

      print(paste("  Matrix:", nrow(P), "x", ncol(P),
                  "| observed links:", sum(P_original)))

      loo_results <- NULL

      for (ri in seq_len(nrow(P))) {
        for (ci in seq_len(ncol(P))) {
          P[ri, ci] <- NA  # withhold this entry

          # C = P with one entry withheld.
          # Equivalent to the self-prediction (/2) branch: when A = P_original,
          # (A + P_with_NA) / 2 reduces to P_original for non-withheld entries
          # and NA for the withheld entry — the same as setting C <- P directly.
          C <- P
          C <- biScale(C, row.center = TRUE, col.center = TRUE,
                       row.scale = FALSE, col.scale = FALSE)
          row_centers <- attr(C, "biScale:row")$center
          col_centers <- attr(C, "biScale:column")$center

          lam0 <- lambda0(C)
          fit  <- softImpute(C, rank.max = 2, lambda = lam0, type = "svd", maxit = 600)
          C_rec      <- softImpute::complete(C, fit)
          C_rec_orig <- C_rec +
            outer(row_centers[rownames(C)], col_centers[colnames(C)], "+")

          loo_results <- rbind(loo_results, data.frame(
            node_to          = rownames(P)[ri],
            node_from        = colnames(P)[ci],
            original_links   = P_original[ri, ci],
            predicted_values = C_rec_orig[ri, ci]
          ))

          P[ri, ci] <- P_original[ri, ci]  # restore entry before next iteration
        }
      }

      combined_results <- rbind(
        combined_results,
        cbind(data.frame(method = m, focal_site = focal_site), loo_results)
      )
    }
  }
  saveRDS(combined_results, results_file)
}

## ---- 4. Classification and evaluation ----
df <- combined_results %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

### ---- Select optimal threshold (max F0.5) ----
# Every entry is a leave-one-out prediction so the full data set can be used
# for threshold evaluation without a separate held-out set.
thresholds <- seq(0, 1, by = 0.1)

df_prepped <- df %>%
  mutate(
    predicted_prob  = sigmoid(predicted_values),
    original_binary = if_else(original_links > 0, 1, 0)
  )

df_thresh <- df_prepped %>%
  tidyr::expand_grid(threshold = thresholds) %>%
  mutate(predicted_bin = if_else(predicted_prob > threshold, 1, 0)) %>%
  group_by(method, focal_site, threshold) %>%
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

df_classified <- df_prepped %>%
  mutate(
    predicted_bin_sigm = if_else(predicted_prob > best_discrete_threshold, 1, 0)
  )

result_summary <- df_classified %>%
  group_by(method, focal_site) %>%
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
  ungroup()

head(result_summary)
summary(result_summary)

## ---- 5. Save combined prediction CSV ----
dir.create("results/predictions", showWarnings = FALSE, recursive = TRUE)

df_out <- df_classified %>%
  group_by(method, focal_site) %>%
  mutate(link_ID = row_number() - 1L) %>%
  ungroup() %>%
  transmute(
    link_ID        = link_ID,
    higher_level   = node_to,                        # pollinator
    lower_level    = node_from,                      # plant
    ground_truth   = as.integer(original_binary),
    probability    = predicted_prob,
    prediction = as.integer(predicted_bin_sigm),
    focal_site,
    method
  )

write_csv(df_out, "results/predictions/serra_marin_loo_prediction_results.csv")
print("Saved: results/predictions/serra_marin_loo_prediction_results.csv")