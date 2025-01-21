# Plot PCA - Figure 2
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

# Define a function to generate scatter plots (correlation)
plot_correlation <- function(data, x_var, y_var, output_file, color_var, line_color, x_label, y_label) {
  ggscatter(data, x = x_var, y = y_var,
                    add = "reg.line",  # Add regression line
                    add.params = list(color = line_color, fill = "lightgray"),
                    cor.coef = TRUE, color = color_var,
                    cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n")) +
    scale_color_npg() +  
    theme(legend.position = "right", axis.text = element_text(size = 6)) +
    xlab(x_label) + 
    ylab(y_label)
}
