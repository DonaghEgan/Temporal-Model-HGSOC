library(tidyverse)
source("src/file_reader.R")

# Load data with error handling
processed_prot_path <- "data/processed/processed_prot.Rds"
meta_data_path <- "data/meta/meta_data.Rds"
gene_meta_path <- "data/meta/gene_meta.Rds"

data <- safe_read_rds(processed_prot_path)
meta <- safe_read_rds(meta_data_path)
gene_meta <- safe_read_rds(gene_meta_path)

# Merge protein data with gene metadata
# Assumes 'data' has rownames as protein IDs and 'gene_meta' has 'Accession' column
data_gene <- data %>% data.frame() %>%
  rownames_to_column(var = "ProteinID") %>%
  left_join(gene_meta, by = c("ProteinID" = "Accession")) %>%
  select(-ProteinID)

# Check for missing gene symbols and warn if any
missing_genes <- data_gene %>%
  filter(is.na(`Gene Symbol`)) %>%
  pull(`Gene Symbol`) %>%
  unique()

if (length(missing_genes) > 0) {
  warning(paste("There are", length(missing_genes), "proteins without gene symbols. They will be removed."))
  data_gene <- data_gene %>%
    filter(!is.na(`Gene Symbol`))
}

# Average duplicated gene symbols by taking the mean
data_gene_avg <- aggregate(. ~ `Gene Symbol`, data = data_gene, FUN = mean) %>%
  column_to_rownames("Gene Symbol")

# Convert the averaged data to a matrix
data_matrix <- as.matrix(data_gene_avg)

# Perform Consensus TME Analysis
data_tme <- ConsensusTME::consensusTMEAnalysis(bulkExp = data_matrix,
  cancerType = "OV",       # Specify the cancer type
  statMethod = "ssgsea"    # Specify the statistical method
)
