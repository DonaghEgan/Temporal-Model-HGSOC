# Load required libraries
library(phenopath)
library(ggplot2)
library(ggpubr)

# ====================================================
# This script performs pseudotime analysis, using phenopath, on  
# the protein expression profiles of multi-site HGSOC tumors. 
# ====================================================

# File paths for input and output
input_data_path <- "data/processed/processed_prot.Rds"
plex_info_path <- "data/meta/meta_data.Rds"
gene_meta_path <- "data/meta/gene_meta.Rds"

# Check if input files exist
if (!file.exists(input_data_path) || 
    !file.exists(plex_info_path) || 
    !file.exists(gene_meta_path)) {
  stop("One or more input files are missing. Please ensure the files are placed in the 'outputs' directory.")
}

# Load data from RDS files
data <- readRDS(input_data_path)
plex_info <- readRDS(plex_info_path)
gene_meta <- readRDS(gene_meta_path)

# Format sample names
plex_info$sample_no <- gsub("S", "P", plex_info$sample_no)

# Run Phenopath
fit <- phenopath(t(data), plex_info$sample_no, elbo_tol = 1e-6, thin = 40)

# Add to meta data
plex_info$pseudotime <- trajectory(fit)

# ====================================================
# Plot results for pseudotime scores
# ====================================================

pdf("analysis/outputs/pseudotime_mets.pdf", height = 3, width = 3)
ggplot(plex_info, aes(x = reorder(metastatic, pseudotime), y = pseudotime)) +
  geom_boxplot() +
  geom_line(aes(group = sample_no), color = "gray", linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.75)) +
  theme_classic() +
  labs(title = "Pseudotime: Metastatic", x = "", y = "Pseudotime") +
  stat_compare_means()
dev.off()

pdf("analysis/outputs/pseudotime_sitegrouped.pdf", height = 4, width = 4)
ggplot(plex_info, aes(x = reorder(site_grouped, pseudotime, mean), y = pseudotime)) +
  geom_jitter(width = 0.1, height = 0, size = 2, aes(color = site_grouped), alpha = 0.7) +  # Add jitter to avoid overplotting
  stat_summary(fun = mean, geom = "crossbar", size = 0.1, color = "red") +  # Add median points
  theme_classic() +
  scale_color_npg() +
  stat_compare_means(method = "anova") +
  labs(title = "Pseudotime: Tumor Site", x = "", y = "Pseudotime") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90))  # Remove legend for color aesthetic
dev.off()

# ====================================================
# Score Pathways
# ====================================================

# Setting up pathways
msigdbr_df <- msigdbr(species = "human", category = "H")
pathways = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

# Clean up row names and unused columns
data_score <- merge(data, gene_meta, by.x = 0, by.y = "Accession")
data_score <- data_score %>%
  dplyr::select(-Row.names) %>%    
  relocate(`Gene Symbol`)   

# Aggregate data by 'Gene Symbol', computing the mean for duplicate genes
data_score <- data_score %>%
  group_by(`Gene Symbol`) %>%
  summarise(across(everything(), mean)) %>%
  na.omit() %>%
  column_to_rownames("Gene Symbol")

# Run GSVA 
data_score <- gsvaParam(as.matrix(data_score), pathways)
data_score <- gsva(data_score)
data_score <- data.frame(t(data_score))

# add pseudotime variable
data_score$pseudotime <- plex_info$pseudotime

# ====================================================
# Determine correlation and plot
# ====================================================

# Select the variable of interest
var_of_interest <- "pseudotime"

# Calculate correlation and p-value
results <- data.frame(variable = character(), correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

for (col in colnames(data_score)) {
  if (col != var_of_interest) {
    test <- cor.test(data_score[[var_of_interest]], data_score[[col]])
    results <- rbind(results, data.frame(variable = col, correlation = test$estimate, p_value = test$p.value))
  }
}

# FDR correction
results$FDR <- p.adjust(results$p_value, method = 'BH', n = nrow(results))

# Top 9 most signifcant
results_top <- results %>% arrange(FDR) %>% dplyr::slice(1:9)
results_top$variable <- gsub("HALLMARK_", "", results_top$variable)
results_top$variable <- gsub("_", " ", results_top$variable)

# Plot
pdf("analysis/outputs/correlations_variable_pseudotime.pdf", height = 2.3, width = 4)
ggplot(results_top, aes(x = reorder(variable, correlation), y = correlation)) +
  geom_col(fill = "skyblue", color = "black", width = 0.4) +
  coord_flip() +
  labs(title = paste("Correlation of", var_of_interest, "with Other Variables"),
       x = "",
       y = "Correlation") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))
dev.off()
