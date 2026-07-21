# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(plotly)
# d3heatmap is not on CRAN yet, but can be found here: https://github.com/talgalili/d3heatmap
# library(d3heatmap)
# 
# # Load data
# data <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/multivariate.csv", header = T, sep = ";")
# colnames(data) <- gsub("\\.", " ", colnames(data))
# 
# # Select a few country
# data <- data %>%
#   filter(Country %in% c("France", "Sweden", "Italy", "Spain", "England", "Portugal", "Greece", "Peru", "Chile", "Brazil", "Argentina", "Bolivia", "Venezuela", "Australia", "New Zealand", "Fiji", "China", "India", "Thailand", "Afghanistan", "Bangladesh", "United States of America", "Canada", "Burundi", "Angola", "Kenya", "Togo")) %>%
#   arrange(Country) %>%
#   mutate(Country = factor(Country, Country))
# 
# # Matrix format
# mat <- data
# rownames(mat) <- mat[,1]
# mat <- mat %>% dplyr::select(-Country, -Group, -Continent)
# mat <- as.matrix(mat)

library(dplyr)
library(tidyr)

# Step 1: Assign numeric codes to link classes
# df_encoded <- df_categorized_with_avg %>%
#   mutate(link_code = as.integer(factor(link_category)))  # 1 = feasible, 2 = forbidden, etc.
# 
# # Step 2: Create wide matrix: node_from as rows, node_to as columns
# link_matrix <- df_encoded %>%
#   select(node_from, node_to, link_code) %>%
#   pivot_wider(names_from = node_to, values_from = link_code) %>%
#   column_to_rownames("node_from")  # matrix rows = plants
# 
# # Step 3: Convert to matrix (with NAs)
# mat <- as.matrix(link_matrix)
# 
# # Optional: Keep key for plotting legend
# link_key <- df_encoded %>%
#   distinct(link_category, link_code) %>%
#   arrange(link_code)

df_interactive_distinct <- distinct(df_categorized_with_avg)

unique(df_interactive_distinct$link_category)

link_matrix <- df_interactive_distinct %>% 
  select(node_from, node_to, link_category) %>% 
  pivot_wider(names_from = node_to, values_from = link_category) %>%
  column_to_rownames("node_from")  # matrix rows = plants

# check if the same interaction gets different category assignments

# Count distinct classes per interaction
conflicts <- df_interactive_distinct %>%
  distinct(node_from, node_to, link_category) %>%
  count(node_from, node_to) %>%
  filter(n > 1)  # Keep only interactions with >1 class

# Optional: View the actual conflicting classes
conflict_classes <- df_interactive_distinct %>%
  semi_join(conflicts, by = c("node_from", "node_to")) %>%
  select(node_from, node_to, link_category) %>%
  distinct() %>%
  arrange(node_from, node_to)

# if it's okay, clean it
df_clean <- df_interactive_distinct %>%
  count(node_from, node_to, link_category) %>%
  group_by(node_from, node_to) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()


# Assign numeric codes
df_clean <- df_clean %>%
  mutate(link_code = as.integer(factor(link_category)))

# Pivot using numeric codes
link_matrix_numeric <- df_clean %>%
  select(node_from, node_to, link_code) %>%
  pivot_wider(names_from = node_to, values_from = link_code) %>%
  column_to_rownames("node_from") %>%
  as.matrix()

# mat <- as.matrix(df_clean)

# heatmap

library(heatmaply)

heatmaply(
  link_matrix_numeric,
  dendrogram = "none",
  xlab = "", ylab = "",
  main = "",
  scale = "none",
  margins = c(60, 100, 40, 20),
  grid_color = "white",
  grid_width = 0.00001,
  titleX = FALSE,
  hide_colorbar = FALSE,
  branches_lwd = 0.1,
  label_names = c("Plant", "Pollinator", "Link class"),
  fontsize_row = 5, fontsize_col = 5,
  labCol = colnames(link_matrix_numeric),
  labRow = rownames(link_matrix_numeric),
  heatmap_layers = theme(axis.line = element_blank()),
  col = viridis::viridis(length(unique(na.omit(as.vector(link_matrix_numeric)))))
)


