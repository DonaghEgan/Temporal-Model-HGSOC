plot_pca <- function(PCAvalues, percentVar, cat_var, colors) {
  ggplot(PCAvalues, aes_string(x = "PC1", y = "PC2", colour = cat_var)) +
    geom_point(size = 4) +
    scale_color_manual(values = colors) +
    theme_bw() +
    theme(legend.text = element_text(size = 8)) +
    labs(
      x = paste0("PC1: ", round(percentVar[1] * 100), "% variance"),
      y = paste0("PC2: ", round(percentVar[2] * 100), "% variance")
    )
}
