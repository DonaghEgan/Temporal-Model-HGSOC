# Load required libraries
library(openxlsx)
library(readr)
library(GSVA)
library(phenopath)
library(ggpubr)
library(ggsci)
library(dplyr)
library(car)
library(ggthemes)

# ====================================================
# Functions
# ====================================================

ttest_two <- function(b1, b2, se_b1, se_b2, n) {
  # Calculate t-statistic
  t_stat <- (b1 - b2) / sqrt((se_b1^2 / n) + (se_b2^2 / n))
  
  # Degrees of freedom
  deg_free <- n - 2
  
  # Calculate one-sided p-values
  p_val_1 <- pt(t_stat, df = deg_free)         # b1 less than b2
  p_val_2 <- 1 - pt(t_stat, df = deg_free)     # b1 greater than b2
  
  # Calculate two-sided p-value
  p_val <- 2 * min(p_val_1, p_val_2)
  
  return(as.numeric(p_val))
}

create_custom_plot <- function(data, x_var, y_var) {
  ggplot(data, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "#2980B9", fill = "grey", alpha = 0.2) +
    labs(
      title = "",
      x = "CD8 Memory:Exhausted",
      y = y_var
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 10)
    ) +
    stat_cor(
      method = "pearson", 
      label.x = min(data[[x_var]], na.rm = TRUE), 
      label.y = max(data[[y_var]], na.rm = TRUE), 
      size = 3.5
    )
}

get_stars <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("")
  }
}

# ===================================================


# Define file paths
input_path <- "/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/inputs/"
output_path <- "/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/"

# Load sample data
sample_number <- read.xlsx(paste0(input_path, "Clincal data IHC patients.xlsx"), sheet = 1) %>%
  distinct(X2, .keep_all = TRUE)

# Load IHC data and clean 'Sample' column
IHC_data <- read.xlsx(paste0(input_path, "file_show (3).xlsx"), sheet = 2) %>%
  mutate(Sample = gsub("\\s+", "", Sample))

# Merge and clean up duplicate rows
IHC_sample <- merge(IHC_data, sample_number, by.x = "Patient.ID", by.y = "X2") %>%
  distinct()

# Group samples based on 'Sample' column
IHC_sample <- IHC_sample %>%
  mutate(Sample_grouped = case_when(
    Sample %in% c("LymphNode", "Spleen") ~ "Lymphoid Mets",
    Sample == "Ovary" ~ "Ovary",
    TRUE ~ "Non Lymphoid Met"
  )) %>%
  filter(!Sample %in% c("Vagina", "Diaphragm")) %>%  # Remove certain samples
  mutate(sample_base = case_when(
    Sample %in% c("LymphNode", "Spleen") ~ "LYMPHOID",
    Sample == "Ovary" ~ "OVARY",
    Sample == "Omentum" ~ "OMENTUM",
    Sample == "Bowel" ~ "BOWEL",
    Sample == "Peritoneum" ~ "PERITONEUM",
    TRUE ~ NA_character_  # Handle missing cases if necessary
  ))

# ===================================================
# Plot correlations
# ===================================================

# Log-transform CD4 and CD8 columns
IHC_sample <- IHC_sample %>%
  mutate(CD4_iTIL_log = log2(CD4_iTIL),
         CD4_sTIL_log = log2(CD4_sTIL),
         CD8_iTIL_log = log2(CD8_iTIL),
         CD8_sTIL_log = log2(CD8_sTIL))

pdf("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/cor_cd4.pdf", height = 3, width = 4)
ggscatter(IHC_sample, x = "CD4_iTIL_log", y = "CD4_sTIL_log",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"),
          cor.coef = TRUE, color = "sample_base",
          cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n")) +
  scale_color_npg() +  theme(legend.position = "right", axis.text = element_text(size = 6)) +
  xlab("Log(CD4 iTIL)") + ylab("Log(CD4 sTIL)")
dev.off()

pdf("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/cor_cd8.pdf", height = 3, width = 4)
ggscatter(IHC_sample, x = "CD8_iTIL_log", y = "CD8_sTIL_log",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "red", fill = "lightgray"),
          cor.coef = TRUE, color = "sample_base",
          cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n")) +
  scale_color_npg() +  theme(legend.position = "right", axis.text = element_text(size = 6)) +
  xlab("Log(CD8 iTIL)") + ylab("Log(CD8 sTIL)")
dev.off()

# ===================================================
# Calculate tan and theta
# ===================================================

IHC_sample$R_CD8 <- sqrt(IHC_sample$CD8_iTIL**2 + IHC_sample$CD8_sTIL**2)
IHC_sample$tan_CD8 <- atan(IHC_sample$CD8_sTIL / IHC_sample$CD8_iTIL)

IHC_sample$R_CD4 <- sqrt(IHC_sample$CD4_iTIL**2 + IHC_sample$CD4_sTIL**2)
IHC_sample$tan_CD4 <- atan(IHC_sample$CD4_sTIL / IHC_sample$CD4_iTIL)

pdf("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/cd8_tan_R.pdf", height = 3, width = 4)
ggplot(IHC_sample, aes(x=R_CD8, y=tan_CD8, color = sample_base)) + geom_point() + theme_bw() + 
  ylab("CD8+ T cells spatial Distribution ()") + xlab("CD8+ T cell quantity (R)") + 
  scale_color_npg() + geom_vline(xintercept = 1000, linetype = "dashed", color = "#ef8a62") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#ef8a62")# Add vertical line
dev.off()

pdf("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/cd4_tan_R.pdf", height = 3, width = 4)
ggplot(IHC_sample, aes(x=R_CD4, y=tan_CD4, color = sample_base)) + geom_point() + theme_bw() + 
  ylab("CD4+ T cells spatial Distribution ()") + xlab("CD4+ T cell quantity (R)")+
  scale_color_npg() 
dev.off()

# ===================================================
# Correlation spatial and pseudotime
# ===================================================

## tme and meta
data_tme <- readRDS("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/data_tme.Rds")
plex_info <- readRDS("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/plex_info.Rds")

## proteomics
data_gs <- readRDS("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/data_gs.Rds")

## sade scored
sade_data <- readRDS("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/sade_data.Rds")

## pheno data
fit <- readRDS("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/inputs/phenopath_fit.Rds")

# Transpose data for further processing
data_gs_t <- as.data.frame(t(data_gs))
data_tme_t <- as.data.frame(t(data_tme))

# Create 'sample_mrn' identifiers for both IHC_sample and plex_info
IHC_sample <- IHC_sample %>%
  mutate(sample_mrn = paste(as.factor(MRN), sample_base, sep = "_"))

plex_info <- plex_info %>%
  mutate(sample_mrn = paste(as.factor(MRN), site_grouped, sep = "_"))

# Add relevant data from transposed data frames
plex_info <- plex_info %>%
  mutate(
    tregs = data_tme_t$T_regulatory_cells,
    cd4 = data_tme_t$T_cells_CD4,
    cd8 = data_tme_t$T_cells_CD8,
    ex_cd8 = sade_data$ex_cd8_tcell,
    mem = sade_data$mem_tcell,
    cytox = sade_data$cytox_lympho,
    ex_cc = sade_data$lymph_ex_cellcycle,
    pseudotime = trajectory(fit))

# Filter plex_info based on matching sample_mrn in IHC_sample
plex_sub <- plex_info %>%
  filter(sample_mrn %in% IHC_sample$sample_mrn)

# Merge plex_sub and IHC_sample based on 'sample_mrn' and remove duplicates
merged_data <- plex_sub %>%
  inner_join(IHC_sample, by = "sample_mrn", multiple = "first")

# Fit linear models for CD8 and CD4 on pseudotime
t_lm_cd8 <- lm(pseudotime ~ scale(R_CD8) + scale(tan_CD8), data = merged_data)
t_lm_cd4 <- lm(pseudotime ~ scale(R_CD4) + scale(tan_CD4), data = merged_data)

# Print summaries
summary(t_lm_cd8)
summary(t_lm_cd4)

# Extract coefficients and standard errors from both models
coef_vals <- data.frame(
  gene = c("CD4_R", "CD4_tan", "CD8_R", "CD8_tan"),
  coef_val = c(
    summary(t_lm_cd4)$coefficients[2, 1],  # CD4_R coefficient
    summary(t_lm_cd4)$coefficients[3, 1],  # CD4_tan coefficient
    summary(t_lm_cd8)$coefficients[2, 1],  # CD8_R coefficient
    summary(t_lm_cd8)$coefficients[3, 1]   # CD8_tan coefficient
  ),
  se_val = c(
    summary(t_lm_cd4)$coefficients[2, 2],  # CD4_R standard error
    summary(t_lm_cd4)$coefficients[3, 2],  # CD4_tan standard error
    summary(t_lm_cd8)$coefficients[2, 2],  # CD8_R standard error
    summary(t_lm_cd8)$coefficients[3, 2]   # CD8_tan standard error
  ),
  ci_lower = c(
    confint(t_lm_cd4)[2, 1],  # CD4_R lower 95% CI
    confint(t_lm_cd4)[3, 1],  # CD4_tan lower 95% CI
    confint(t_lm_cd8)[2, 1],  # CD8_R lower 95% CI
    confint(t_lm_cd8)[3, 1]   # CD8_tan lower 95% CI
  ),
  ci_upper = c(
    confint(t_lm_cd4)[2, 2],  # CD4_R upper 95% CI
    confint(t_lm_cd4)[3, 2],  # CD4_tan upper 95% CI
    confint(t_lm_cd8)[2, 2],  # CD8_R upper 95% CI
    confint(t_lm_cd8)[3, 2]   # CD8_tan upper 95% CI
  ),
  p_val = c(
    summary(t_lm_cd4)$coefficients[2, 4],  # CD4_R p-value
    summary(t_lm_cd4)$coefficients[3, 4],  # CD4_tan p-value
    summary(t_lm_cd8)$coefficients[2, 4],  # CD8_R p-value
    summary(t_lm_cd8)$coefficients[3, 4]
))


# Add stars based on p-values
coef_vals$stars <- cut(
  coef_vals$p_val,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "")
)

pdf("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/coefficients_R_tan_pseudotime.pdf", height = 2.5, width = 2)
ggplot(coef_vals, aes(x = gene, y = coef_val)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, color = gene), width = 0) +
  geom_point(aes(colour = gene), size = 3) +
  geom_point(shape = 21, size = 3, fill = "white") +
  labs(x = "", y = "FCI") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        axis.title.y = element_text(size = 12)) +
  ylab("Pseudotime") +
  xlab("") +
  scale_color_manual(values = c("#1661ec", '#1661ec', '#f23030', '#f23030')) +
  # Add stars above the points
  geom_text(aes(label = stars, y = coef_val + 0.1), color = "black", size = 7)
dev.off()

# ===================================================
# CD4 / Tregs
# ===================================================

# Fit linear model
t_lm_cd4 <- lm(scale(R_CD4) ~ scale(cd4), data = merged_data)
t_lm_tregs <- lm(scale(R_CD4) ~ scale(tregs), data = merged_data)

summary(t_lm_cd4)
summary(t_lm_tregs)

# Display summary of the model
summary(t_lm_tregs)
summary(t_lm_cd4)

# Extract coefficients and standard errors from both models

coef_vals <- data.frame(
  gene = c("Tregs", "CD4"),
  coef_val = c(
    summary(t_lm_tregs)$coefficients[2, 1],  # Coefficient for Tregs
    summary(t_lm_cd4)$coefficients[2, 1]     # Coefficient for CD4
  ),
  se_val = c(
    summary(t_lm_tregs)$coefficients[2, 2],  # Standard error for Tregs
    summary(t_lm_cd4)$coefficients[2, 2]     # Standard error for CD4
  ),
  ci_lower = c(
    confint(t_lm_tregs)[2, 1],  # Lower 95% CI for Tregs
    confint(t_lm_cd4)[2, 1]     # Lower 95% CI for CD4
  ),
  ci_upper = c(
    confint(t_lm_tregs)[2, 2],  # Upper 95% CI for Tregs
    confint(t_lm_cd4)[2, 2]     # Upper 95% CI for CD4
  )
)

# Create plot for coefficient comparison with error bars
pdf("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/coefficients_tregs_cd4_R.pdf", height = 2.5, width = 2)
ggplot(coef_vals, aes(x = gene, y = coef_val)) + 
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, color = gene), width = 0) +
  geom_point(aes(colour = gene), size = 3) +
  geom_point(shape = 1,size = 3,colour = "black") + 
  labs(x = "", y = "FCI") +  
  theme_classic() +
  ylim(-.2, 1.3) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") + ylab("Pseudotime") + xlab("") +
  scale_color_manual(values = c("#f29f05", "#1661ec"))
dev.off()

p_val_tregs <- ttest_two(coef_vals$coef_val[1], coef_vals$coef_val[2], coef_vals$se_val[1], coef_vals$se_val[2], n = nrow(merged_data))

# ===================================================
# CD8s
# ===================================================

t_lm_mem  <- lm(scale(pseudotime) ~ scale(mem) + scale(cytox) + scale(ex_cd8) + scale(tan_CD8), data = merged_data)
summary(t_lm_mem)

# Extract summary statistics from the model
coef_summary <- summary(t_lm_mem)$coefficients[-1, ]  # Remove intercept

# Create a dataframe with coefficients, standard errors, and p-values
coef_vals <- data.frame(
  term = rownames(coef_summary),
  coef_val = coef_summary[, "Estimate"],
  se_val = coef_summary[, "Std. Error"],
  p_val = coef_summary[, "Pr(>|t|)"]
)

# Calculate 95% confidence intervals
conf_int <- confint(t_lm_mem)[-1, ]  # Remove intercept
coef_vals$ci_lower <- conf_int[, 1]
coef_vals$ci_upper <- conf_int[, 2]

# Apply the function to assign stars
coef_vals$stars <- sapply(coef_vals$p_val, get_stars)

# Optional: Rename terms for better readability
coef_vals$term <- dplyr::recode(coef_vals$term,
                         "scale(mem)" = "Memory CD8⁺",
                         "scale(cytox)" = "Cytotoxic CD8⁺",
                         "scale(ex_cd8)" = "Exhausted CD8⁺",
                         "scale(tan_CD8)" = "Tangent CD8⁺")

# Define colors for each term
term_colors <- c("Memory CD8⁺" = "#22BABB",
                 "Cytotoxic CD8⁺" = "#F39C12",
                 "Exhausted CD8⁺" = "#3498DB",
                 "Tangent CD8⁺" = "#9B59B6")  # Add more colors if needed

# Create plot for coefficient comparison with error bars
pdf("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/coefficients_ex_mem_cd8_tan.pdf", height = 2.8, width = 2)
ggplot(coef_vals, aes(x = term, y = coef_val)) + 
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, color = term), width = 0) +
  geom_point(aes(colour = term), size = 3) +
  geom_point(shape = 21, size = 3, fill = "white") +  # Add a border to points
  geom_text(aes(label = stars), vjust = -1.5, color = "black", size = 5) +  # Add stars above points
  labs(x = "", y = "Pseudotime") +  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        axis.title.y = element_text(size = 12)) + 
  scale_color_manual(values = term_colors) + 
  ylim(min(coef_vals$ci_lower) - 0.1, max(coef_vals$ci_upper) + 0.1)
dev.off()
