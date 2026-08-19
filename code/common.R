
# Shared helpers for the operational pipeline. Trimmed to only what
# code/serra_marin_svd_prediction_loo.R actually calls (sigmoid,
# build_interaction_matrix) — everything else (Canary/emln-specific helpers,
# plotting functions, the old bootstrap-imputation routine) lived here for
# the legacy scripts in code/legacy/ and isn't needed by the current
# pipeline. See code/legacy/ if those are ever needed again.

# functions -------------

# transforming raw predictions for binary evaluation
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

build_interaction_matrix <- function(data, layers_to_filter) {

  # unify layer format
  layers <- layers_to_filter # default, they are character
  if (is.numeric(layers_to_filter)) layers <- paste0("layer_", layers_to_filter)

  # Step 1: Filter rows based on specified layers
  filtered_data <- subset(data, layer_from %in% layers)

  # Step 2: Aggregate weights for identical species pairs

  aggregated_data <- filtered_data %>%
    group_by(node_from, node_to) %>%
    summarise(weight = sum(weight), .groups = 'drop')

  # Step 3: Create the matrix with specific row and column species
  species_from <- unique(aggregated_data$node_from)  # Columns
  species_to <- unique(aggregated_data$node_to)      # Rows

  # Initialize an empty matrix
  interaction_matrix <- matrix(0, nrow = length(species_to), ncol = length(species_from),
                               dimnames = list(species_to, species_from))

  # Populate the matrix with aggregated weights
  for (i in 1:nrow(aggregated_data)) {
    row <- aggregated_data$node_to[i]    # Rows represent 'node_to' species
    col <- aggregated_data$node_from[i]  # Columns represent 'node_from' species
    interaction_matrix[row, col] <- aggregated_data$weight[i]
  }

  return(interaction_matrix)
}
