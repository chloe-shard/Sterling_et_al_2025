################################ IVY GAP analysis ##############################
################################################################################

# Required libraries
library(readr)
library(dplyr)

#load Ivy Gap data
IVY_GAP_Matrix <- read_csv("file location to IVY_GAP_Matrix.csv")

#specifiy genes of interest
Genesofinterest <- c("DHODH", "PRMT5", "WDR77")

#subset matrix for genes of interest
Genesofinterest_IVY_GAP <- subset(IVY_GAP_Matrix, subset = Region %in% Genesofinterest)
Genesofinterest_IVY_GAP <- data.frame(Genesofinterest_IVY_GAP, row.names = 1)

# Ensure the object is a data frame
Genesofinterest_IVY_GAP <- as.data.frame(Genesofinterest_IVY_GAP)

# Order the columns alphabetically and select them
Genesofinterest_IVY_GAP <- Genesofinterest_IVY_GAP %>% dplyr::select(order(colnames(Genesofinterest_IVY_GAP)))

print(Genesofinterest_IVY_GAP)

## Calculate the average expression per region

# Extract the matrix
expression_matrix <- Genesofinterest_IVY_GAP

# Get the column names and extract region names
region_names <- gsub("\\.\\.\\..*", "", colnames(expression_matrix))

# Create a data frame to hold the averaged results
averaged_expression <- data.frame(row.names = rownames(expression_matrix))

# Calculate mean expression for each region
for (region in unique(region_names)) {
  # Subset columns for the current region
  region_columns <- expression_matrix[, region_names == region, drop = FALSE]
  
  # Calculate the row means
  averaged_expression[[region]] <- rowMeans(region_columns, na.rm = TRUE)
}

# Display the averaged expression matrix
print(averaged_expression)

##Plot a heatmap
# Required library
library(pheatmap)

# Define the desired gene and region order
gene_order <- c("DHODH", "PRMT5", "WDR77")
region_order <- c("Leading.Edge", "Infiltrating.Tumor", "Cellular.Tumor", 
                  "Hyperplastic.blood.vessels", "Microvascular.proliferation", 
                  "Perinecrotic.zone", "Pseudopalisading.cells.around.necrosis")

# Reorder the rows and columns of the averaged expression matrix
# Ensure that only the specified genes and regions are used
heatmap_data <- averaged_expression[gene_order, region_order, drop = FALSE]

# Generate the heatmap
pheatmap(heatmap_data,
         cluster_rows = FALSE,     # Do not cluster rows
         cluster_cols = FALSE,     # Do not cluster columns
         color = colorRampPalette(c("blue", "white", "red"))(50), # Colour scale
         fontsize_row = 10,        # Font size for rows (genes)
         fontsize_col = 10)        # Font size for columns (regions)

## z score transform the data across regions

# Perform z-score transformation across the regions (columns)
zscore_transformed <- t(apply(heatmap_data, 1, scale))
colnames(zscore_transformed) <- colnames(heatmap_data)
rownames(zscore_transformed) <- rownames(heatmap_data)

# Generate the heatmap
pheatmap(zscore_transformed,
         cluster_rows = FALSE,     # Do not cluster rows
         cluster_cols = FALSE,     # Do not cluster columns
         color = colorRampPalette(c("blue", "yellow", "red"))(50),
         #color = colorRampPalette(c("blue", "white", "red"))(50),# Colour scale
         fontsize_row = 10,        # Font size for rows (genes)
         fontsize_col = 10)        # Font size for columns (regions)

# Perform Anova analysis to identify if genes are sig upreg in a region

# Required libraries
library(dplyr)
library(multcomp)

# Extract region names
region_names <- gsub("\\.\\.\\..*", "", colnames(expression_matrix))

# Initialize a list to store ANOVA p-values
anova_results_list <- list()

# Initialize a list to store Tukey HSD results
tukey_results_list <- list()

# Perform ANOVA and Tukey HSD test for each gene
for (gene in rownames(expression_matrix)) {
  
  # Extract gene expression data as a numeric vector
  gene_expression <- as.numeric(expression_matrix[gene, ])
  
  # Create a dataframe for the current gene
  df <- data.frame(Expression = gene_expression, Region = factor(region_names))
  
  # Perform ANOVA
  aov_res <- aov(Expression ~ Region, data = df)
  p_value <- summary(aov_res)[[1]][["Pr(>F)"]][1]  # Extract p-value
  
  # Store p-value
  anova_results_list[[gene]] <- p_value
  
  # Perform Tukey HSD post hoc test if ANOVA is significant (p < 0.05)
  if (p_value < 0.05) {
    tukey_res <- TukeyHSD(aov_res, "Region")
    
    # Convert Tukey results to a dataframe
    tukey_df <- as.data.frame(tukey_res$Region)
    tukey_df$Gene <- gene  # Assign correct gene name
    tukey_df$Comparison <- rownames(tukey_df)  # Add comparison column
    rownames(tukey_df) <- NULL  # Reset row names
    
    # Store in the list
    tukey_results_list[[gene]] <- tukey_df
  }
}

# Convert ANOVA results to a dataframe
anova_df <- data.frame(Gene = names(anova_results_list),
                       P_value = unlist(anova_results_list))

# Adjust p-values using Benjamini-Hochberg correction
anova_df$Adjusted_P_value <- p.adjust(anova_df$P_value, method = "BH")

# Filter significant genes (FDR < 0.05)
significant_genes <- anova_df %>% filter(Adjusted_P_value < 0.05)

# Convert Tukey results list into a single dataframe
if (length(tukey_results_list) > 0) {
  tukey_results_df <- do.call(rbind, tukey_results_list)
} else {
  tukey_results_df <- data.frame()  # Empty dataframe if no significant results
}


