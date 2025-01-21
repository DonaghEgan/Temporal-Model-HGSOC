library(ggpubr)
source("src/plot_functions.R")
source("src/stat_test_functions.R")

# ===================================================
# Load Data
# ===================================================

ihc_data <- read.xlsx("data/processed/Supplemental_Tables.xlsx", sheet = 4)

# ===================================================
# Plot correlations
# ===================================================

# Log-transform CD4 and CD8 columns
IHC_sample <- ihc_data %>%
  mutate(CD4_iTIL_log = log2(CD4_iTIL),
         CD4_sTIL_log = log2(CD4_sTIL),
         CD8_iTIL_log = log2(CD8_iTIL),
         CD8_sTIL_log = log2(CD8_sTIL))

p1 <- plot_correlation(IHC_sample, "CD4_iTIL_log", "CD4_sTIL_log", 
                       color_var = "sample_base", line_color ="blue",
                       x_label = "Log(CD4 iTIL)",y_label = "Log(CD4 sTIL)")

p2 <- plot_correlation(IHC_sample, "CD8_iTIL_log", "CD8_sTIL_log", 
                       color_var = "sample_base", line_color ="red",
                       x_label = "Log(CD8 iTIL)",y_label = "Log(CD8 sTIL)")

pdf("analysis/outputs/cor_iTILs_sTILs.pdf", onefile = T)
p1
p2
dev.off()

# ===================================================
# Calculate tan and R
# ===================================================

IHC_sample$R_CD8 <- sqrt(IHC_sample$CD8_iTIL**2 + IHC_sample$CD8_sTIL**2)
IHC_sample$tan_CD8 <- atan(IHC_sample$CD8_sTIL / IHC_sample$CD8_iTIL)

IHC_sample$R_CD4 <- sqrt(IHC_sample$CD4_iTIL**2 + IHC_sample$CD4_sTIL**2)
IHC_sample$tan_CD4 <- atan(IHC_sample$CD4_sTIL / IHC_sample$CD4_iTIL)

p3 <- ggplot(IHC_sample, aes(x=R_CD8, y=tan_CD8, color = sample_base)) + geom_point() + theme_bw() + 
  ylab("CD8+ T cells spatial Distribution ()") + xlab("CD8+ T cell quantity (R)") + 
  scale_color_npg()

p4 <- ggplot(IHC_sample, aes(x=R_CD4, y=tan_CD4, color = sample_base)) + geom_point() + theme_bw() + 
  ylab("CD4+ T cells spatial Distribution ()") + xlab("CD4+ T cell quantity (R)")+
  scale_color_npg() 

pdf("analysis/outputs/polar_coordinates.pdf", onefile = T)
p3
p4
dev.off()

# ===================================================
# Correlation spatial and pseudotime
# ===================================================

## deconvolution data and meta
data_tme <- readRDS("data/processed/data_tme.Rds")
meta <- readRDS("data/meta/meta_data.Rds")

## proteomics with gene identifiers
data_gs <- readRDS("data/processed/data_gs.Rds")

## sade signatures scored in proteomics
sade_data <- readRDS("data/processed/sade_data.Rds")

## pseudotime estimates
fit <- readRDS("data/processed/phenopath_fit.Rds")

# Transpose data for further processing
data_gs_t <- as.data.frame(t(data_gs))
data_tme_t <- as.data.frame(t(data_tme))

# Create 'sample_mrn' identifiers for both IHC_sample and plex_info
IHC_sample <- IHC_sample %>%
  mutate(sample_mrn = paste(as.factor(MRN), sample_base, sep = "_"))

meta <- meta %>%
  mutate(sample_mrn = paste(as.factor(MRN), site_grouped, sep = "_"))

# Add relevant data from transposed data frames
meta <- meta %>%
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
plex_sub <- meta %>%
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

pdf("analysis/outputs/coefficients_R_tan_pseudotime.pdf", height = 2.5, width = 2)
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
pdf("analysis/outputs/coefficients_ex_mem_cd8_tan.pdf", height = 2.8, width = 2)
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

