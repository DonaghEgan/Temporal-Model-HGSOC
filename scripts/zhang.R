# Load required libraries
library(ppcor)
library(DESeq2)
library(tidyverse)
library(corrplot)
library(GSVA)
library(readxl)

# Load count data and metadata
load("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/inputs/EGAD00010001515.RData")
counts_zhang <- EGAD00010001515_raw_count_tb
counts_meta <- EGAD00010001515_metadata

# Format counts table by removing unnecessary columns and adjusting column names
counts_zhang <- counts_zhang %>%
  column_to_rownames("Name") %>%
  dplyr::select(-CodeClass, -Accession, -Accession_nv) # Remove unnecessary columns

colnames(counts_zhang) <- sub("_[^_]*$", "", colnames(counts_zhang)) # Simplify column names

# Merge with metadata for correct sample IDs
counts_zhang <- counts_zhang %>%
  t() %>%
  merge(counts_meta[, c("ns_sample_id", "sample_id")], by.x = 0, by.y = "ns_sample_id") %>%
  column_to_rownames("sample_id") %>%
  dplyr::select(-Row.names) # Remove redundant row names column

# DESeq2 normalization
dds <- DESeqDataSetFromMatrix(countData = t(counts_zhang), colData = counts_meta, design = ~ 1)
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)

# Extract rlog-transformed values
rld <- rlog(dds, blind = FALSE)
normalized_data <- assay(rld)

# Load metadata with cluster identities
zhang_meta <- read_xlsx("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/inputs/zhang_meta.xlsx", sheet = 2)

# Subset genes significantly up in infiltrated and excluded samples based on limma results
up_inf <- result_dep %>%
  filter(logFC > 0, P.Value < 0.05) %>%
  pull(`Gene Symbol`)

ex_sig <- result_dep %>%
  filter(logFC < 0, P.Value < 0.05) %>%
  pull(`Gene Symbol`)

# Create gene sets for GSVA scoring
gene_sets <- list(inf_sig = up_inf, ex_sig = ex_sig)

# GSVA scoring on normalized data
score_zhang <- GSVA::zscoreParam(normalized_data, geneSets = gene_sets)
score_zhang <- GSVA::gsva(score_zhang)

# Merge GSVA scores with zhang_meta
zhang_scores_merged <- merge(t(score_zhang), zhang_meta, by.x = 0, by.y = "sample_id")

# Remove NAs from merged data
zhang_scores_merged <- zhang_scores_merged %>%
  filter(!is.na(til_cluster))

# Define comparisons for statistical tests
my_comp <- list(c("ES-TIL", "N-TIL"), c("ES-TIL", "S-TIL"))

# Create color scheme for plots
til_colors <- c("red", "grey", "blue")

# PDF for infiltrated genes - Infiltrated Signature
pdf("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/zhang_cluster_inf_genes.pdf")
ggplot(zhang_scores_merged, aes(x = reorder(til_cluster, inf_sig, FUN = median), y = inf_sig, fill = til_cluster)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comp) +
  scale_fill_manual(values = til_colors) +
  theme_bw() +
  xlab("") +
  ylab("Infiltrated Signature")
dev.off()

# PDF for excluded genes - Excluded Signature
pdf("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/zhang_cluster_ex_genes.pdf")
ggplot(zhang_scores_merged, aes(x = reorder(til_cluster, -ex_sig, FUN = median), y = ex_sig, fill = til_cluster)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comp) +
  scale_fill_manual(values = til_colors) +
  theme_bw() +
  xlab("") +
  ylab("Excluded Signature")
dev.off()

################################################################################
# Correlation with Immune Characteristics

# Load metadata
zhang_meta <- read_xlsx("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/inputs/zhang_meta.xlsx", sheet = 8)

# Merge score_zhang data with metadata and convert patient_id to numeric
zhang_s_m <- merge(t(score_zhang), zhang_meta, by.x = 0, by.y = "sample_id") %>%
  mutate(patient_id = as.numeric(patient_id))

# Select relevant columns for correlation analysis
ex_sig <- zhang_s_m[, 2:8]

# Define variables for partial correlation
target_variable <- "inf_sig"
control_variable <- "patient_id"
other_variables <- setdiff(names(ex_sig), c(target_variable, control_variable))
variables <- setdiff(names(ex_sig), control_variable)

## CORRELATION MATRIX

# Initialize matrices for partial correlations and p-values
partial_corr_matrix <- matrix(NA, nrow = length(variables), ncol = length(variables),
                              dimnames = list(variables, variables))
partial_p_matrix <- matrix(NA, nrow = length(variables), ncol = length(variables),
                           dimnames = list(variables, variables))

# Calculate pairwise partial Spearman correlations and p-values
for (i in seq_along(variables)) {
  for (j in seq_along(variables)) {
    if (i != j) {
      partial_result <- pcor.test(ex_sig[[variables[i]]], ex_sig[[variables[j]]], ex_sig[[control_variable]], method = "spearman")
      partial_corr_matrix[i, j] <- partial_result$estimate
      partial_p_matrix[i, j] <- partial_result$p.value
    }
  }
}

# Cap p-values below 0.01
partial_p_matrix[partial_p_matrix < 0.01] <- 0.01

# Save the correlation plot to a PDF
pdf("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/zhang_inf_signature_immune_metrics.pdf")

corrplot(partial_corr_matrix, 
         method = "color", col = rev(COL2('RdBu', 200)), 
         type = "lower", 
         tl.col = "black", tl.srt = 45,   # Text label color and rotation
         p.mat = partial_p_matrix, sig.level = 0.01, pch.cex = 1.5, insig = "p-value", 
         diag = FALSE)  
         
dev.off()

################################################################################

# Filter and select the relevant columns from ex_sig
ex_sig <- zhang_s_m[, c(2,3,4,11:14)]

# Variables to be used in partial correlation calculation
control_variable <- "patient_id"
variables <- setdiff(names(ex_sig), control_variable)

# Initialize matrices for partial correlation coefficients and p-values
partial_corr_matrix <- matrix(NA, nrow = length(variables), ncol = length(variables),
                              dimnames = list(variables, variables))
partial_p_matrix <- matrix(NA, nrow = length(variables), ncol = length(variables),
                           dimnames = list(variables, variables))

# Calculate pairwise partial Spearman correlations and p-values
for (i in seq_along(variables)) {
  for (j in seq_along(variables)) {
    if (i != j) {
      partial_corr <- pcor.test(ex_sig[[variables[i]]], ex_sig[[variables[j]]], ex_sig[[control_variable]], method = "spearman")
      partial_corr_matrix[i, j] <- partial_corr$estimate
      partial_p_matrix[i, j] <- partial_corr$p.value
    }
  }
}

partial_p_matrix[partial_p_matrix < 0.01] <- 0.01
corrplot(partial_corr_matrix, method="color", col = rev(COL2('RdBu', 200)), 
         type="lower", 
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = partial_p_matrix, sig.level = 0.0, pch.cex = 1.5, insig = "p-value", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)
