# Load necessary libraries
library(tidyverse)
library(readxl)
library(lmtest)
library(ggsci)
library(msigdbr)
library(fgsea)
source("src/stat_test_functions.R")

# ====================================================
# This script validates relationship between Keratin expression
# and T cell infiltration and further validations in Sanchez et al
# ====================================================

# Define file paths relative to the current working directory
meta_data_path <- "data/processed/41588_2020_630_MOESM3_ESM.xlsx"
norm_counts_path <- "data/processed/TreatmentNaive_log2exp_data.txt"

# Read metadata: Skip the first 2 rows, and convert relevant columns to numeric
meta_data <- read_xlsx(
  path = meta_data_path,    # Use relative paths
  sheet = 4,                # Specify the Excel sheet
  skip = 2                  # Skip the first 2 rows
) %>%
  mutate(
    # Convert specific columns to numeric
    across(c(CD8, CD4, FoxP3, CELL_COUNTS), as.numeric),
    
    # Compute fractions
    CD8_frac   = CD8   / CELL_COUNTS,
    CD4_frac   = CD4   / CELL_COUNTS,
    FoxP3_frac = FoxP3 / CELL_COUNTS,
    T_frac     = (CD4 + CD8 + FoxP3) / CELL_COUNTS
  )

# Read normalized counts and set "Hugo_Symbol" as row names
norm_counts <- read.delim(
  file = norm_counts_path,
  header = TRUE,         # Ensure the first row contains column names
  stringsAsFactors = FALSE  # Avoid converting strings to factors
) %>%
  column_to_rownames("Hugo_Symbol")  # Move "Hugo_Symbol" column to row names

# Remove rows with NA values
complete_data <- na.omit(meta_data)

# Compute average metrics per WELL
case_data <- complete_data %>%
  group_by(WELL) %>%
  summarise(
    CD8_frac   = mean(CD8_frac,   na.rm = TRUE),
    CD4_frac   = mean(CD4_frac,   na.rm = TRUE),
    T_frac     = mean(T_frac,     na.rm = TRUE),
    FoxP3_frac = mean(FoxP3_frac, na.rm = TRUE),
    Total_counts = mean(CELL_COUNTS, na.rm = TRUE)
  )

# Merge case_data with transposed norm_counts by matching "WELL" to row names
case_norm <- merge(
  case_data,
  t(norm_counts),
  by.x = "WELL",
  by.y = 0
)

# Compute mean KRT expression:
# Reverse log2 by raising 2 to the power of expression, average them, then take log2 again.
case_norm$mean_krt <- log2((case_norm$KRT1^2 + case_norm$KRT2^2 +
                              case_norm$KRT9^2 + case_norm$KRT10^2) / 4)

# Convert fractions to percentages
case_norm <- case_norm %>%
  mutate(
    T_frac     = T_frac * 100,
    CD8_frac   = CD8_frac * 100,
    CD4_frac   = CD4_frac * 100,
    FoxP3_frac = FoxP3_frac * 100
  )

# Fit linear models
ifi44l_lm <- lm(T_frac ~ IFI44L, data = case_norm)
summary(ifi44l_lm)

krt_lm <- lm(T_frac ~ mean_krt, data = case_norm)
summary(krt_lm)

# Get coefficients from models
coef_1 <- summary(krt_lm)$coefficients[2,1]
coef_2 <- summary(ifi44l_lm)$coefficients[2,1]
se_b1 = summary(krt_lm)$coefficients[2,2]
se_b2 = summary(ifi44l_lm)$coefficients[2,2]
n = nrow(case_norm)

# Compare using T test
p_val <- ttest_two(coef_1, coef_2, se_b1, se_b2, n)

# ====================================================
# Plot results for coefficients
# ====================================================

# Define coefficients and standard errors
coef_vals <- c(coef_1, coef_2)   
se_vals   <- c(se_b1, se_b2)     

# Create a data frame for coefficients with confidence intervals
coef_df <- data.frame(
  gene    = c("KRT", "IFI44L"),
  coef_val = coef_vals,
  se_val   = se_vals
)

# Calculate confidence intervals using confint()
ci_krt     <- confint(krt_lm)
ci_ifi44l  <- confint(ifi44l_lm)

# Add confidence intervals to the data frame (2nd row is the slope)
coef_df$ci_lower <- c(ci_krt[2, 1], ci_ifi44l[2, 1])  
coef_df$ci_upper <- c(ci_krt[2, 2], ci_ifi44l[2, 2])


# Generate the plot
pdf("analysis/outputs/coefficients_krt_ifi44l.pdf", height = 2.5, width = 1.8)
ggplot(coef_df, aes(x = gene, y = coef_val, color = gene)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0) +
  geom_point(size = 3) +
  geom_point(shape = 1, size = 3, color = "black") +
  labs(x = "", y = "T cell Infiltration") +
  theme_classic() +
  theme(
    axis.text.x      = element_text(angle = 90),
    legend.position  = "none"
  ) +
  scale_color_manual(values = c("#56B4E9", "#66cc99"))
dev.off()

# =====================================================
# Co-expression Analysis - KRT
# =====================================================

# Load and process normalized counts
norm_counts <- read.delim("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/inputs/TreatmentNaive_log2exp_data.txt") %>%
  column_to_rownames("Hugo_Symbol") %>%
  data.frame() %>%
  t() %>%
  as.data.frame()

# Compute mean KRT expression and remove individual KRT columns
norm_counts <- norm_counts %>%
  mutate(
    mean_krt = log2((KRT1^2 + KRT2^2 + KRT9^2 + KRT10^2) / 4)
  ) %>%
  dplyr::select(-KRT1, -KRT2, -KRT9, -KRT10)

# Calculate correlations with mean_krt for all other genes
correlations <- data.frame(cor = sapply(norm_counts, function(x) cor(norm_counts$mean_krt, x))) %>%
  rownames_to_column("gene") %>%
  filter(gene != "mean_krt") %>%
  mutate(rank = rank(-cor, ties.method = "first"))

# Set up pathways using the msigdbr package
msigdbr_df <- msigdbr(species = "human", category = "C5")
pathways <- lapply(split(msigdbr_df$gene_symbol, msigdbr_df$gs_name), unique)

# Convert correlations to named vector for fgsea analysis
genes <- correlations %>%
  arrange(cor) %>%
  dplyr::select(gene, cor) %>%
  deframe()

# Perform fgsea analysis
fgseaRes <- fgseaSimple(pathways = pathways, stats = genes, nperm = 100000)

# Identify TCR-related genes
tcr_genes <- pathways[["GOCC_T_CELL_RECEPTOR_COMPLEX"]]
top_5_percent_threshold <- ceiling(0.05 * nrow(correlations))

# Assign labels for plotting TCR genes
correlations <- correlations %>%
  mutate(
    label = ifelse(gene %in% tcr_genes & rank < top_5_percent_threshold, gene, "NONE")
  )

# plot correlations
pdf("analysis/outputs/coex_krt_genes.pdf")
ggplot(correlations) +
  geom_point(aes(x = rank, y = cor, color = label), size = 0.5) +
  geom_point(data = subset(correlations, !label == "NONE"),
             aes(x = rank, y = cor, color = label)) +
  geom_vline(xintercept = top_5_percent_threshold, linetype = "dashed", color = "grey") +  # Add vertical line
  labs(x = "Rank", y = "Correlation") +  # Add axis labels and title
  theme_bw() + scale_color_manual(values = c("#56B4E9","#000000", "#56B4E9", "#56B4E9"))
dev.off()

# plot enrichment
pdf("analysis/outputs/tcr_complex_krt_coexp.pdf")
plotEnrichment(pathways[["GOCC_T_CELL_RECEPTOR_COMPLEX"]],
               genes) + labs(title = "GOCC_T_CELL_RECEPTOR_COMPLEX")
dev.off()

# =====================================================
# Correlation KRT w/different cell types
# =====================================================

# Test different cell types with linear models
tests <- list(
  FOXP3 = lm(scale(FoxP3_frac) ~ mean_krt, data = case_norm),
  CD4   = lm(scale(CD4_frac) ~ mean_krt, data = case_norm),
  CD8   = lm(scale(CD8_frac) ~ mean_krt, data = case_norm)
)

# Extract coefficients, standard errors, and confidence intervals
coef_df <- do.call(rbind, lapply(names(tests), function(cell_type) {
  model <- tests[[cell_type]]
  coef_val <- summary(model)$coefficients[2, 1]
  se_val <- summary(model)$coefficients[2, 2]
  ci <- confint(model)[2, ]
  data.frame(
    CELL_TYPE = cell_type,
    COEF_VAL = coef_val,
    SE_VAL = se_val,
    CI_LOWER = ci[1],
    CI_UPPER = ci[2]
  )
}))

# Optionally, calculate p-values for specific comparisons
n <- nrow(case_norm)
p_val_cd4_foxp3 <- ttest_two(coef_df$COEF_VAL[coef_df$CELL_TYPE == "FOXP3"],
                             coef_df$COEF_VAL[coef_df$CELL_TYPE == "CD4"],
                             coef_df$SE_VAL[coef_df$CELL_TYPE == "FOXP3"],
                             coef_df$SE_VAL[coef_df$CELL_TYPE == "CD4"], n)

p_val_cd8_foxp3 <- ttest_two(coef_df$COEF_VAL[coef_df$CELL_TYPE == "FOXP3"],
                             coef_df$COEF_VAL[coef_df$CELL_TYPE == "CD8"],
                             coef_df$SE_VAL[coef_df$CELL_TYPE == "FOXP3"],
                             coef_df$SE_VAL[coef_df$CELL_TYPE == "CD8"], n)

# Print the coefficient table for verification
print(coef_df)

# Plot the results
pdf("analysis/outputs/coefficients_cell_type_krt.pdf", height = 2, width = 2.5)
ggplot(coef_df, aes(x = CELL_TYPE, y = COEF_VAL, color = CELL_TYPE)) +
  geom_errorbar(aes(ymin = CI_LOWER, ymax = CI_UPPER), width = 0) +
  geom_point(size = 3) +
  geom_point(shape = 1, size = 3, color = "black") +
  coord_flip() +
  labs(x = "", y = "KRT Expression") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#56B4E9", "#66CC99", "#FF3300"))
dev.off()                     