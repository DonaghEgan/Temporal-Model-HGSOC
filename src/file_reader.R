# Function to safely load RDS files with error handling
safe_read_rds <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("File does not exist:", file_path))
  }
  readRDS(file_path)
}
