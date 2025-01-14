# Load necessary libraries
library(ConsensusTME)
library(ggpubr)
library(gridExtra)
library(readxl)
library(openxlsx)
library(ggsci)
library(ppcor)
library(ggrepel)
library(RColorBrewer)
library(tidyverse)
library(GSVA)
library(matrixStats)
library(infer)

# Generate color palettes
set3_colors <- brewer.pal(6, "Set3")
dark2_colors <- brewer.pal(6, "Dark2")
pastel1_colors <- brewer.pal(6, "Pastel1")
colors <- c(set3_colors, dark2_colors, pastel1_colors)

# Read proteomics data
data <- readRDS("outputs/imputed_data.Rds")
plex_info <- readRDS("outputs/plex_info.Rds")
gene_meta <- readRDS("outputs/gene_meta.Rds")

# Process plex info
plex_info <- plex_info %>%
  mutate(
    site_lymph = case_when(
      site_grouped %in% "LYMPHOID" ~ "Lymphoid Mets",
      site_grouped %in% "OVARY" ~ "Ovary",
      TRUE ~ "Non Lymphoid Met"
    )
  )

# Merge and aggregate gene data
data_gs <- merge(data, gene_meta, by.x = 0, by.y = "Accession") %>%
  select(-Row.names) %>%
  group_by(`Gene Symbol`) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("Gene Symbol")

data_tme <- readRDS("/outputs/data_tme.Rds")

# Filter out non-cell types
df <- data_tme[!rownames(data_tme) %in% c("Immune_Score"), ]

# Function to calculate coefficient of variation (CV)
calc_cv <- function(row) {
  sd(row) / mean(row + 1) # Avoid division by zero with small adjustment
}

cv_values <- data.frame(cv = abs(apply(df, 1, calc_cv)), cell_type = rownames(df))
cv_values$rank <- rank(dplyr::desc(cv_values$cv))

# Plot coefficient of variation
pdf("outputs/coef_variation_cell_types.pdf", height = 3, width = 3)
ggplot(cv_values, aes(x = rank, y = cv, label = cell_type)) +
  geom_point() +
  geom_label_repel(size = 2) +
  labs(y = "Coefficient of Variation (CV)", x = "Rank") +
  theme_bw()
dev.off()

# Permutation testing
melt_df <- reshape2::melt(df)
colnames(melt_df)[1] <- "cell_type"

# Define parameters
sample_size <- nrow(df)
perm_reps <- 100000

many.perm <- melt_df %>%
  rep_sample_n(size = sample_size, replace = FALSE, reps = perm_reps) %>%
  mutate(perm_treatment = sample(cell_type, size = n(), replace = FALSE)) %>%
  group_by(replicate, perm_treatment)

# Calculate CV for permutations
many.perm.means <- many.perm %>%
  group_by(replicate) %>%
  summarise(cv_cell = sd(value) / mean(value + 1), .groups = "drop")

many.perm.diffs <- many.perm.means %>%
  summarise(diff_cv = diff(cv_cell))

# Plot permuted sampling distribution
pdf("outputs/test.pdf")
ggplot(many.perm.diffs, aes(x = diff_cv)) +
  geom_histogram(bins = 32, color = "white") +
  labs(title = "Permuted Sampling Distribution", x = "Difference in CV") +
  theme_bw()
dev.off()

# Observed differences
observed <- cv_values %>%
  filter(cell_type %in% c("T_regulatory_cells", "T_cells_CD8", "T_cells_CD4", "T_cells_gamma_delta")) %>%
  pivot_wider(names_from = cell_type, values_from = cv)

observed_diffs <- observed %>%
  summarise(
    Tregs_CD4 = T_regulatory_cells - T_cells_CD4,
    Tregs_gamma = T_regulatory_cells - T_cells_gamma_delta,
    Tregs_CD8 = T_regulatory_cells - T_cells_CD8,
    CD4_CD8 = T_cells_CD8 - T_cells_CD4
  ) %>%
  pivot_longer(cols = everything(), names_to = "comparison", values_to = "value")

# P-value calculation function
calculate_p_value <- function(observed_diff) {
  (sum(abs(many.perm.diffs$diff_cv) >= abs(observed_diff)) + 1) / (nrow(many.perm.diffs) + 1)
}

observed_diffs <- observed_diffs %>%
  mutate(p_values = sapply(value, calculate_p_value))

# Null distribution with observed differences
mean_diff <- mean(many.perm.diffs$diff_cv)
sd_diff <- sd(many.perm.diffs$diff_cv)
ci_lower <- mean_diff - 2 * sd_diff
ci_upper <- mean_diff + 2 * sd_diff

# Plot histogram with mean and 95% CI
pdf("outputs/different_cv_null.pdf", height = 2.5, width = 3)
ggplot(many.perm.diffs, aes(x = diff_cv)) +
  geom_histogram(bins = 32, fill = "#bebebe", color = "white") +
  geom_vline(xintercept = ci_lower, color = "red", linetype = "dashed", size = 0.7) +
  geom_vline(xintercept = ci_upper, color = "red", linetype = "dashed", size = 0.7) +
  geom_vline(xintercept = 0.12511825, color = colors[7], linetype = "dashed", size = 0.7) +
  geom_vline(xintercept = 0.10682240, color = colors[6], linetype = "dashed", size = 0.7) +
  geom_vline(xintercept = 0.11290375, color = "#0073c2", linetype = "dashed", size = 0.7) +
  
  labs(
    title = "Null Distribtion",
    x = "Difference in CV"
  ) +
  theme_bw()
dev.off()
