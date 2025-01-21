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
