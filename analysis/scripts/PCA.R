library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(msigdbr)
library(GSVA)
source("src/plot_functions.R")

# Generate 6 colors from each palette
dark2_colors <- pal_nejm()(8)
pastel1_colors <- pal_jco()(10)
pastel2_colors <- pal_lancet()(9)
npg_colors <- pal_npg()(9)

# Combine the colors from all three palettes
colors <- c(dark2_colors, pastel1_colors, pastel2_colors)

# Read proteomics 
data <- readRDS("data/processed/processed_prot.Rds")
meta <- readRDS("data/meta/meta_data.Rds")
meta$sample_no <- gsub("S", "P", meta$sample_no)

# Set NA values to zero
data[is.na(data)] <- 0

# Running PCA 
pca <- prcomp(t(data), scale. = T, center = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

# Add meta to PCA
PCAvalues$metastatic <- meta$metastatic 
PCAvalues$site <- meta$sample_base 
PCAvalues$sample <- as.factor(meta$sample_no)
PCAvalues$site_grouped <- meta$site_grouped

# plot
p1 <- plot_pca(PCAvalues, percentVar, colors = npg_colors, cat_var = "site_grouped")
p2 <- plot_pca(PCAvalues, percentVar, colors = c(pastel1_colors, pastel2_colors, npg_colors), cat_var = "sample")
p3 <- plot_pca(PCAvalues, percentVar, colors = dark2_colors, cat_var = "metastatic")

# Open a PDF device
pdf("analysis/outputs/pca_files.pdf", onefile = T)
p1
p2
p3
dev.off()
