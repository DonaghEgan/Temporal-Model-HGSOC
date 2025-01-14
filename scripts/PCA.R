library(ggplot2)
library(viridis)
library(ggsci)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(msigdbr)
library(GSVA)

# Generate 6 colors from each palette
dark2_colors <- pal_nejm()(8)
pastel1_colors <- pal_jco()(10)
pastel2_colors <- pal_lancet()(9)

# Combine the colors from all three palettes
colors <- c(dark2_colors, pastel1_colors, pastel2_colors)

# read proteomics 
data <- readRDS("outputs/imputed_data.Rds")
plex_info <- readRDS("outputs/plex_info.Rds")
plex_info$sample_no <- gsub("S", "P", plex_info$sample_no)

# set NA values to zero
data[is.na(data)] <- 0

# Running PCA 
pca <- prcomp(t(data), scale. = T, center = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

# add meta
PCAvalues$metastatic <- plex_info$metastatic 
PCAvalues$site <- plex_info$sample_base 
PCAvalues$sample <- as.factor(plex_info$sample_no)
PCAvalues$site_grouped <- plex_info$site_grouped

# plot
p1 <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = site_grouped)) +
  geom_point(size = 4) + scale_color_viridis(discrete = T) + theme_bw() + theme(legend.text = element_text(size =8)) +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))

p2 <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = site)) +
  geom_point(size = 4) + scale_color_manual(values = colors) + theme_bw() + theme(legend.text = element_text(size =8)) +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))

p3 <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = sample)) +
  geom_point(size = 4) + scale_color_manual(values = colors) + theme_bw() + theme(legend.text = element_text(size =8)) +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))

p4 <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = metastatic)) +
  geom_point(size = 4) + scale_color_manual(values = colors) + theme_bw() + theme(legend.text  = element_text(size =8)) +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))

# Open a PDF device
pdf("outputs/protypia_pca.pdf", onefile = T)
p1
p2
p3
p4
dev.off()