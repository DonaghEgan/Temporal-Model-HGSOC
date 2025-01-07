# Load necessary libraries
library(entropy)
library(ggpubr)
library(readr)
library(readxl)
library(tidyverse)
library(reshape2)
library(msigdbr)
library(dplyr)
library(glmnet)

# ====================================================
# This script calculates the entropy for CD4+ and CD8+ TILs and 
# determines what pathways explain the variation in entropy using 
# an elastic net
# ====================================================

# ====================================================
# Functions
# ====================================================

# Function to calculate Shannon entropy for a vector
calculate_entropy <- function(x) {
  
  # Normalize the counts to probabilities
  freq <- x / sum(x)
  
  entropy_val <- entropy::entropy(freq, method = "ML")
  
  return(entropy_val)
}

# ===================================================

# Read clinical data
clinical_data <- read_xlsx("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/inputs/Clincal data IHC patients.xlsx",
                           sheet = 1, trim_ws = TRUE)
# Read in IHC data
ihc_data <- read_xlsx("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/inputs/Clincal data IHC patients.xlsx",
                      sheet = 3)

# Merge data and remove duplicates
ihc_sample <- ihc_data %>%
  left_join(clinical_data[, c(1, 2)], by = c("Sample" = "...2"),multiple = "first") 

# Determine number of each phenotype - CD4 
cd4_counts <- ihc_sample %>%
  dplyr::select(Type, MRN, Ovary, Omentum) %>%
  pivot_longer(cols = c(Ovary, Omentum), names_to = "site", values_to = "Count") %>%
  mutate(MRN_site = paste(MRN, site, sep = "_")) %>%
  group_by(Type, MRN_site) %>%
  summarise(Total_Count = sum(Count)) %>%
  pivot_wider(names_from = MRN_site, values_from = Total_Count) %>%
  column_to_rownames("Type")

# Determine number of each phenotype - CD8 
cd8_counts <- ihc_sample %>%
  dplyr::select(Type, MRN, `Ovary CD8`, `Omentum CD8`) %>%
  pivot_longer(cols = c(`Ovary CD8`, `Omentum CD8`), names_to = "site", values_to = "Count") %>%
  mutate(MRN_site = paste(MRN, site, sep = "_")) %>%
  group_by(Type, MRN_site) %>%
  summarise(Total_Count = sum(Count)) %>%
  pivot_wider(names_from = MRN_site, values_from = Total_Count) %>%
  column_to_rownames("Type")


# Calculate entropy
entropy_cd4 <- data.frame(entropy_cd4 = apply(cd4_counts, 2, calculate_entropy),
                          sample = toupper(colnames(cd4_counts))) %>% melt()

entropy_cd8 <- data.frame(entropy_cd8 = apply(cd8_counts, 2, calculate_entropy),
                          sample = gsub(" CD8", "",toupper(colnames(cd8_counts)))) %>% melt()

# Combine entropy
entrop_comp <- rbind(entropy_cd4, entropy_cd8)

## plot
pdf("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/compare_entropy.pdf", height = 2.5, width = 2.5)
ggplot(entrop_comp, aes(x=reorder(variable, value), y = value, fill = variable)) + geom_boxplot(alpha = 0.75) + 
  stat_compare_means() + theme_classic() +
  scale_fill_manual(values = c("#67a9cf", "#8a62ef")) + xlab("") + ylab("Shannon Entropy") +
  scale_x_discrete(labels = c("CD8","CD4")) +
  theme(legend.position = "none")
dev.off()

# ===================================================
# Score pathways
# ===================================================

# Load the data
protein_expression_data <- readRDS("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/data_gs.Rds")
plex_info <- readRDS("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/plex_info.Rds")

# Transpose the data and set appropriate row names
protein_expression_trans <- data.frame(t(protein_expression_data))
rownames(protein_expression_trans) <- plex_info$`Protypia ID`

# Load Pathways
msigdbr_df <- msigdbr(species = "human", category = "H")
pathways = lapply(split(msigdbr_df$gene_symbol, msigdbr_df$gs_name), unique)

# Score pathways 
pathways_score <- ssgseaParam(as.matrix(protein_expression_data), pathways, minSize = 10)
pathways_score <- gsva(pathways_score)
pathways_score <- data.frame(t(pathways_score))

# Set row names using Protypia ID from plex_info and merge
rownames(pathways_score) <- plex_info$`Protypia ID`
data_meta <- merge(plex_info, pathways_score, by.x = "Protypia ID", by.y = 0)

# Merge plex_info and protein data, filter for specific samples, and create site_mrn
data_meta <- data_meta %>%
  filter(sample_base %in% c("OMENTUM", "OVARY")) %>%  
  mutate(site_mrn = paste(as.factor(MRN), site_grouped, sep = "_"))

# Filter data_meta for matched sites
matched_sites <- intersect(data_meta$site_mrn,toupper(colnames(cd4_counts)))
data_meta <- data_meta %>%
  filter(site_mrn %in% matched_sites) 
 
# Combine entropy values for CD8 and CD4
entropy_tme <- bind_cols(entropy_cd8, entropy_cd4)

# format
entropy_tme <- entropy_tme %>%
  dplyr::rename(cd8_ent = value...3, cd4_ent = value...6, sample = sample...1) %>%
  dplyr::select(sample, cd8_ent, cd4_ent) 

# Join entropy and meta data
data_meta <- merge(entropy_tme, data_meta, by.x = "sample", by.y = "site_mrn")
  
# Subet pathways and entropy columns
data_entropy <- data_meta[,colnames(data_meta) %in% c("cd8_ent", "cd4_ent", colnames(pathways_score))]

# ===================================================
# Elastic Net - CD8
# ===================================================

# Create independent and target variables
X <- as.matrix(data_entropy[, which(names(data_entropy) %in% colnames(pathways_score))])  # all columns except 'y' (target)
y <- data_entropy$cd8_ent  # target variable

# CV
set.seed(123)  # For reproducibility
lasso_model <- cv.glmnet(X, y, alpha = 1, grouped = F, thresh = 1e-3, maxit = 10e+10)

# select best lambda 
best_lambda <- lasso_model$lambda.min

# Plot cross-validation results
plot(lasso_model)
# Add a vertical line at log(0.01)

# Fit model
best_model <- glmnet(X, y, alpha = 1, lambda =  best_lambda)

# Select coefficients
coefficients <- as.matrix(coef(best_model))
non_zero_coef <- coefficients[coefficients != 0, , drop = FALSE]

# Format
coef_df <- as.data.frame(non_zero_coef[-1, , drop = FALSE])  # Exclude the intercept
coef_df$Variable <- rownames(coef_df)  # Add variable names as a new column
names(coef_df) <- c("Coefficient", "Variable")

# Plot
pdf("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/lasso_cd8_ent.pdf", height = 2.2, width = 5)
ggplot(coef_df, aes(x = reorder(Variable, Coefficient), y = Coefficient)) +
  geom_bar(stat = "identity", fill = "#8a62ef", width = 0.7, alpha = 0.75, color = "black") +  
  coord_flip() + 
  labs(title = "CD8: Non-zero Coef", x = "", y = "Coefficient") +  
  theme_classic() + 
  theme(axis.text.y = element_text(size = 8),  
        axis.text.x = element_text(size = 6.5),
        plot.title = element_text(hjust = 0.5, size = 10)) 
dev.off()

# ===================================================
# Elastic Net - CD4
# ===================================================

# Reset target variable
y <- data_entropy$cd4_ent  # target variable

set.seed(123)  # For reproducibility
lasso_model <- cv.glmnet(X, y, alpha = 1, grouped = F)

# Plot cross-validation results
plot(lasso_model)
best_lambda <- lasso_model$lambda.min

# Fit model
best_model <- glmnet(X, y, alpha = 1, lambda = best_lambda)

# Extract coefficients
coefficients <- as.matrix(coef(best_model))
non_zero_coef <- coefficients[coefficients != 0, , drop = FALSE]

# Formant output
coef_df <- as.data.frame(non_zero_coef[-1, , drop = FALSE])  # Exclude the intercept
coef_df$Variable <- rownames(coef_df)  # Add variable names as a new column
names(coef_df) <- c("Coefficient", "Variable")

pdf("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/lasso_cd4_ent.pdf", height = 2.2, width = 5)
ggplot(coef_df, aes(x = reorder(Variable, Coefficient), y = Coefficient)) +
  geom_bar(stat = "identity", fill = "#67a9cf", width = 0.7, alpha = 0.75, color = "black") +  
  coord_flip() + 
  labs(title = "CD4: Non-zero Coef", x = "", y = "Coefficient") +  
  theme_classic() + 
  theme(axis.text.y = element_text(size = 8),  
        axis.text.x = element_text(size = 6.5),
        plot.title = element_text(hjust = 0.5, size = 10)) 
dev.off()

# ===================================================
# GC correlation
# ===================================================

# Perform correlation test
cor_test <- cor.test(sub_entropy$cd8_ent, sub_entropy$GC)

# Extract correlation results
correlation <- cor_test$estimate
p_value <- cor_test$p.value

pdf("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/cor_cd8_gc.pdf", height = 2, width = 2)
ggplot(sub_entropy, aes(x = GC, y = cd8_ent)) +
  geom_point(alpha = 0.6, size = 2) +  # Scatter plot points
  geom_smooth(method = "lm", color = "#8a62ef", se = TRUE) +  # Add linear regression line
  labs(
    x = "Vitamin D Binding Protein",
    y = "CD8 SDI"
  ) +
  annotate(
    "text",
    x = min(sub_entropy$GC, na.rm = TRUE),
    y = max(sub_entropy$cd8_ent, na.rm = TRUE),
    label = paste0("r = ", round(correlation, 2), "p = ", signif(p_value, 2)),
    hjust = 0,
    size = 3.2
  ) +
  theme_classic()
dev.off()
