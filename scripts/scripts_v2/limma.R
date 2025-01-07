library(limma)
library(data.table)
library(fgsea)
library(msigdbr)
library(tidyverse)
library(RColorBrewer)
library(EnhancedVolcano)
library(ggpubr)
library(readxl)
library(dplyr)

# Define file paths
data_file <- "/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/imputed_data.Rds"
plex_info_file <- "/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/plex_info.Rds"
gene_meta_file <- "/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/gene_meta.Rds"
clinical_data_file <- "/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/inputs/Clincal data IHC patients.xlsx"

# Load data
proteomics_data <- readRDS(data_file)
plex_info <- readRDS(plex_info_file)
gene_meta <- readRDS(gene_meta_file)

# Set proteomics data column names using MRN and sample base from plex_info
colnames(proteomics_data) <- paste(plex_info$MRN, plex_info$sample_base, sep = "_")

# Load clinical data (sample numbers and IHC data)
sample_number <- read_xlsx(clinical_data_file, sheet = 1)
IHC_data <- read_xlsx(clinical_data_file, sheet = 3)

# Merge IHC data with sample numbers and ensure distinct entries
IHC_sample <- merge(IHC_data, sample_number[, c(1, 2)], by.x = "Sample", by.y = "...2") %>%
  distinct()

# Create combined CD4/CD8 scores for ovary and omentum
IHC_sample <- IHC_sample %>%
  mutate(
    CD4_CD8_ovary = (Ovary + `Ovary CD8`) / 2,
    CD4_CD8_Omentum = (`Omentum CD8` + Omentum) / 2
  )

# Classify based on maximum CD4/CD8 score for ovary
ovary_class <- IHC_sample %>%
  group_by(`Patient ID`) %>%
  dplyr::slice(which.max(CD4_CD8_ovary)) %>%
  dplyr::select(`Patient ID`, Type) %>%
  mutate(Samples = "OVARY")

# Classify based on maximum CD4/CD8 score for omentum
omentum_class <- IHC_sample %>%
  group_by(`Patient ID`) %>%
  dplyr::slice(which.max(CD4_CD8_Omentum)) %>%
  dplyr::select(`Patient ID`, Type) %>%
  mutate(Samples = "OMENTUM")

# Combine ovary and omentum classifications
omentum_ovary_class <- bind_rows(omentum_class, ovary_class)

# Remove specific sample not present in IHC data
omentum_ovary_class <- omentum_ovary_class %>%
  filter(`Patient ID` != "3421487")

# Create unique patient identifier by combining Patient ID and sample type (ovary/omentum)
omentum_ovary_class <- omentum_ovary_class %>%
  mutate(Patient.ID = paste(`Patient ID`, Samples, sep = "_"))

# Subset proteomics data based on matching patient IDs
proteomics_data_sub <- proteomics_data[, colnames(proteomics_data) %in% omentum_ovary_class$Patient.ID]

# Extract components from omentum_ovary_class
Samples <- omentum_ovary_class$Samples
p_id <- sub("_.*", "", omentum_ovary_class$Patient.ID) 
ihc <- omentum_ovary_class$Type

# Create design matrix directly and apply renaming in one step
design <- model.matrix(~ 0 + ihc + p_id + Samples)
colnames(design) <- gsub("Samples|ihc", "", colnames(design))

# Define contrasts for linear modeling
contrast <- limma::makeContrasts(
  inf_vs_non_inf = Infiltrated - (Desert + excluded / 2),
  levels = colnames(design)
)

# Subset gene_meta and data_sub to match relevant genes
gene_meta_sub <- gene_meta[!is.na(gene_meta$`Gene Symbol`),]
proteomics_data_sub <- proteomics_data_sub[rownames(proteomics_data_sub) %in% gene_meta_sub$Accession,]

# Fit linear model and apply contrasts
fit <- lmFit(proteomics_data_sub, design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)

# Extract top results and merge with gene_meta for annotations
result_dep <- topTable(fit, n=Inf)
result_dep <- merge(result_dep, gene_meta, by.x = 0, by.y = "Accession")
saveRDS(result_dep, "/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/DEP_infiltrated.Rds")

pdf("/mnt/8TB/Projects/POIAZ/Donagh/protypia_proteomics/outputs/volcano_inflamed_.pdf", height = 4.5, width = 6,  pointsize = 8)
EnhancedVolcano(result_dep,
                lab = result_dep$`Gene Symbol`, title = "Infiltrated vs Non infiltrated",
                x = "logFC",
                y = "P.Value", pCutoff = 0.05, legendPosition = "top",raster = TRUE,pointSize = 1.5,
                FCcutoff = 1, col = c("grey", "grey", "grey", "red"), labSize = 2, 
                xlim = c(-4,4), ylim = c(0,8)) + theme_classic()
dev.off()