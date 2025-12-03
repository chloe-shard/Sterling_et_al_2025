################### Bulk RNAseq analysis of KDM family genes in GSC from GSE119834 ##################

# Load required packages
library(readr)
library(pheatmap)

# load expression matrix (GSE119834)
expr_matrix <- read.csv("file path/GSE119834_fpkm_table.csv", row.names = 1)

# Specify genes of interest
genes_of_interest <- c("KDM1A", "KDM1B", "KDM2A", "KDM2B", "KDM3A", "KDM3B", "JMJD1C",
                       "KDM4A", "KDM4B", "KDM4C", "KDM4D", "KDM5A", "KDM5B", "KDM5C",
                       "KDM5D", "KDM6A", "KDM6B", "KDM7A", "PHF8", "PHF2", "KDM8")

# Filter rows (genes)
filtered_expr <- expr_matrix[rownames(expr_matrix) %in% genes_of_interest, ]

# Remove samples (columns) starting with "NSC"
filtered_expr <- filtered_expr[, !grepl("^NSC", colnames(filtered_expr))]

# Remove samples (columns) starting with "GBM"
filtered_expr <- filtered_expr[, !grepl("^GBM", colnames(filtered_expr))]

# Log2 transform
log_expr <- log2(filtered_expr + 1)

# Heatmap without row scaling
pheatmap(log_expr,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         scale = "none",
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("#010180", "white", "#EE9337"))(50),
         main = "KDM Family Genes in GSC")
