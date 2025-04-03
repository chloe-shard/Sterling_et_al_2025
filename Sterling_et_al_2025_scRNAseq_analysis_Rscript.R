#################### Single cell Analysis #####################
###############################################################

# Required library
library(Seurat)

####################### UMAPs of malignant cells in Neftel dataset showing co-expression ##########################

#load dataset
obj_Neftel <- readRDS("file location to Neftel.rds")

# subset for adult GBM samples
Idents(obj_Neftel) <- "GBMType"
obj_Neftel <- subset(obj_Neftel, idents = c("Adult"))

# subset for malignant cells
Idents(obj_Neftel) <- "CellAssignment"
obj_Neftel_Mal <- subset(obj_Neftel, idents = c("Malignant"))

# Required libraries
library(Seurat)
library(ggplot2)

## for the following code replace WDR77 with PRMT5 to generate an alternative UMAP for the figure

# 1. Extract expression data for the genes of interest 
gene_expr <- FetchData(obj_Neftel_Mal, vars = c("DHODH", "WDR77"))

# 2. Create a new column 'gene_combo' based on expression:
# A cell is considered positive if its count is > 0.
gene_expr$gene_combo <- ifelse(gene_expr$DHODH > 0 & gene_expr$WDR77 > 0, "DHODH+ WDR77+",
                               ifelse(gene_expr$DHODH > 0 & gene_expr$WDR77 == 0, "DHODH+ WDR77-",
                                      ifelse(gene_expr$DHODH == 0 & gene_expr$WDR77 > 0, "DHODH- WDR77+",
                                             "DHODH- WDR77-")))

# 3. Add the new gene_combo information to the Seurat object's metadata
obj_Neftel_Mal <- AddMetaData(obj_Neftel_Mal, metadata = gene_expr$gene_combo, col.name = "gene_combo")

# 4. Define the colors for each gene combination
color_map <- c("DHODH+ WDR77+" = "#F0900B",
               "DHODH- WDR77-" = "#D3D3D3",
               "DHODH- WDR77+" = "#1F618D",
               "DHODH+ WDR77-" = "#85C1E9")


# 5. Generate the UMAP plot using DimPlot, grouping cells by the new gene_combo column
DimPlot(obj_Neftel_Mal, reduction = "umap", group.by = "gene_combo", cols = color_map) +
  ggtitle("UMAP annotated by DHODH and WDR77 expression")

############### Score Neftel cell states in single cell Seurat object

# Required libraries
library(Seurat)
library(reshape2)

# Read Neftel subgroup marker genes

gbm.genes.tab <- read.table("neftel_markers.txt", header = TRUE, row.names = NULL, sep = "\t")
gbm.genes.tab <- gbm.genes.tab[,!(colnames(gbm.genes.tab) %in% c("G1_S", "G2_M"))]
gbm.genes.lst <- melt(gbm.genes.tab, measure.vars = colnames(gbm.genes.tab))
gbm.genes.lst <- gbm.genes.lst[gbm.genes.lst$value != "",]

gbm.genes.uniq <- unique(gbm.genes.lst$value)

table(gbm.genes.lst$variable)

# Calculate the meta module score for each Neftel subgroup

obj <- readRDS("file location to Obj_Seurat.rds")

Idents(obj) <- "Assignment"
obj.mal <- subset(obj, ident = "Glioma")
obj.mal.df <- GetAssayData(object = obj.mal[["RNA"]], slot = "data")

# Calculate aggregate expression per gene and remove genes with 0 expression in tumor cells
genes.aggregate.expr <- sort(rowMeans(obj.mal.df))
genes.aggregate.expr <- genes.aggregate.expr[genes.aggregate.expr > 0]

obj.mal.df <- obj.mal.df[row.names(obj.mal.df) %in% names(genes.aggregate.expr),]

gbm_states <- unique(gbm.genes.lst$variable)

state_scores <- matrix(NA,nrow=ncol(obj.mal.df),ncol=length(gbm_states),dimnames = list(colnames(obj.mal.df), gbm_states))

for(g in gbm_states){
  
  cat(g, "\n")
  
  genes <- gbm.genes.lst[gbm.genes.lst$variable == g,]$value
  
  obj.mal.state <- obj.mal.df[row.names(obj.mal.df) %in% genes,]
  obj.mal.nonstate <- obj.mal.df[!row.names(obj.mal.df) %in% genes,]
  
  # Calculate relative expression for state
  mean_exp_tumor_state <- rowMeans(obj.mal.state)
  obj.mal.state.rel <- obj.mal.state/mean_exp_tumor_state
  mean_rel_tumor_state <- colMeans(obj.mal.state.rel)
  
  # Calculate relative expression for non-state
  mean_exp_tumor_nonstate <- rowMeans(obj.mal.nonstate)
  obj.mal.nonstate.rel <- obj.mal.nonstate/mean_exp_tumor_nonstate
  mean_rel_tumor_nonstate <- colMeans(obj.mal.nonstate.rel)
  
  
  # Calculate state score
  state.score <- mean_rel_tumor_state/mean_rel_tumor_nonstate
  
  state_scores[, g] <- state.score
}

cellstate <- apply(state_scores[,1:6], 1, function(x) ifelse(max(x) > 1, names(which.max(x)), "None"))
obj.mal <- AddMetaData(obj.mal, metadata = cellstate, col.name = "cellstate")
obj <- AddMetaData(obj, metadata = cellstate, col.name = "cellstate")

### Save the object now for reloading
saveRDS(obj, file = "Obj_Seurat.rds")

####################### Dotplot showing expression of genes of interest across Neftel cell states ##########################

# Required libraries
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)

#load datasets and subset for malignant cells from primary GBM

##Ebert
obj_Ebert_Mal <- readRDS("file location to Ebert.rds")
#subset cells assigned to a cellstate
Idents(obj_Ebert_Mal) <- "cellstate"
obj_Ebert_Mal <- subset(obj_Ebert_Mal, idents = c("AC", "MES1.Hypoxia.independent.", "MES2.Hypoxia.dependent.", "NPC1", "NPC2", "OPC"))


##LeBlanc
obj_LeBlanc <- readRDS("file location to LeBlanc.rds")
#subset to remove recurrent samples
Idents(obj_LeBlanc) <- "sample"
obj_LeBlanc <- subset(obj_LeBlanc, idents = c("JK124_reg1_tis_1", "JK124_reg1_tis_2", "JK124_reg2_tis_1", "JK124_reg2_tis_2", "JK125_reg1_tis_1.1", "JK125_reg1_tis_1.2", "JK125_reg2_tis_1", "JK125_reg2_tis_2_r1", "JK126_reg1_tis_1.1", "JK126_reg1_tis_1.2", "JK126_reg2_tis_1", "JK134_reg1_tis_1", "JK134_reg2_tis_1", "JK136_reg1_tis_1", "JK136_reg2_tis_1", "JK136_reg2_tis_2_br", "JK142_reg1_tis_1", "JK142_reg2_tis_1", "JK142_reg2_tis_2.1_br", "JK142_reg2_tis_2.2_br", "JK152_reg1_tis_1", "JK152_reg2_tis_1", "JK153_reg1_tis_1", "JK153_reg2_tis_1", "JK156_reg1_tis_1", "JK156_reg2_tis_1", "JK156_reg2_tis_2_br", "JK163_reg1_tis_1", "JK163_reg2_tis_1"))
# subset malignant cells
Idents(obj_LeBlanc) <- "cell_type"
obj_LeBlanc_Mal <- subset(obj_LeBlanc, idents = c("malignant"))
#subset cells assigned to a cellstate
Idents(obj_LeBlanc_Mal) <- "cellstate"
obj_LeBlanc_Mal <- subset(obj_LeBlanc_Mal, idents = c("AC", "MES1.Hypoxia.independent.", "MES2.Hypoxia.dependent.", "NPC1", "NPC2", "OPC"))

##Neftel
obj_Neftel <- readRDS("file location to Neftel.rds")
Idents(obj_Neftel) <- "GBMType"
obj_Neftel <- subset(obj_Neftel, idents = c("Adult"))
# subset malignant cells
Idents(obj_Neftel) <- "CellAssignment"
obj_Neftel_Mal <- subset(obj_Neftel, idents = c("Malignant"))
#subset cells assigned to a cellstate
Idents(obj_Neftel_Mal) <- "cellstate"
obj_Neftel_Mal <- subset(obj_Neftel_Mal, idents = c("AC", "OPC", "MES2", "MES1","NPC1","NPC2"))

##Abdelfattah
obj_Abdelfattah <- readRDS("file location to Abdelfattah.rds")
# subset to remove recurrent samples
Idents(obj_Abdelfattah) <- "Type"
obj_Abdelfattah <- subset(obj_Abdelfattah, idents = c("GBM"))
# subset malignant cells
Idents(obj_Abdelfattah) <- "Assignment"
obj_Abdelfattah_Mal <- subset(obj_Abdelfattah, idents = c("Glioma"))
#subset cells assigned to a cellstate
Idents(obj_Abdelfattah_Mal) <- "cellstate"
obj_Abdelfattah_Mal <- subset(obj_Abdelfattah_Mal, idents = c("AC", "MES1_HypoxIndep", "MES2_HypoxDep", "NPC1", "NPC2", "OPC"))

#Chen
obj_Chen <- readRDS("file location to Chen.rds")
# subset to remove recurrent samples
Idents(obj_Chen) <- "orig.ident"
obj_Chen <- subset(obj_Chen, idents = c("PDC001", "PJ052", "PJ053", "PW016-703", "PW017-703", "PW032-710", "PW035-710", "PW039-705"))
# subset malignant cells
Idents(obj_Chen) <- "cellkb"
obj_Chen_Mal <- subset(obj_Chen, idents = c("Tumor"))
#subset cells assigned to a cellstate
Idents(obj_Chen_Mal) <- "cellstate"
obj_Chen_Mal <- subset(obj_Chen_Mal, idents = c("AC", "MES1_HypoxIndep", "MES2_HypoxDep", "NPC1", "NPC2", "OPC"))


# Couturier ####
obj_Couturier <- readRDS("file location to Couturier.rds")
# subset to remove recurrent samples
Idents(obj_Couturier) <- "orig.ident"
obj_Couturier <- subset(obj_Couturier, idents = c("BT333", "BT338_1of2", "BT338_2of2", "BT346", "BT363_1of2", "BT363_2of2", "BT364_1of2", "BT364_2of2", "BT368", "BT389", "BT390", "BT397_1of2", "BT397_2of2", "BT400", "BT402", "BT407", "BT409"))
# subset malignant cells
Idents(obj_Couturier) <- "cellkb"
obj_Couturier_Mal <- subset(obj_Couturier, idents = c("Tumor"))
#subset cells assigned to a cellstate
Idents(obj_Couturier_Mal) <- "cellstate"
obj_Couturier_Mal <- subset(obj_Couturier_Mal, idents = c("AC", "MES1_HypoxIndep", "MES2_HypoxDep", "NPC1", "NPC2", "OPC"))

# Dotplot

# Required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

genes_of_interest <- c("DHODH", "PRMT5", "WDR77")

# Map dataset-specific cell state labels to standardized labels
normalize_cellstate <- function(cellstate, dataset_name) {
  mapping <- list(
    "Neftel" = c("MES1" = "MES1", "MES2" = "MES2"),
    "Ebert" = c("MES1.Hypoxia.independent." = "MES1", "MES2.Hypoxia.dependent." = "MES2"),
    "Adelfattah" = c("MES1_HypoxIndep" = "MES1", "MES2_HypoxDep" = "MES2"),
    "Couturier" = c("MES1_HypoxIndep" = "MES1", "MES2_HypoxDep" = "MES2"),
    "Chen" = c("MES1_HypoxIndep" = "MES1", "MES2_HypoxDep" = "MES2"),
    "LeBlanc" = c("MES1.Hypoxia.independent." = "MES1", "MES2.Hypoxia.dependent." = "MES2")
  )
  
  # Ensure the mapping works for all datasets
  standardized <- mapping[[dataset_name]][cellstate]
  return(ifelse(is.na(standardized), cellstate, standardized))
}

# Function to calculate average expression (mean), z-score, and percentage expression per dataset
process_dataset <- function(seurat_obj, dataset_name) {
  # Normalize cellstate labels
  seurat_obj@meta.data$cellstate <- sapply(seurat_obj@meta.data$cellstate, normalize_cellstate, dataset_name)
  
  # Filter for genes that are present in the dataset
  genes_present <- intersect(genes_of_interest, rownames(seurat_obj))
  
  # Handle cases where no genes are present
  if (length(genes_present) == 0) {
    stop(paste("No genes of interest found in the dataset:", dataset_name))
  }
  
  # Extract gene expression matrix
  assay_data <- GetAssayData(seurat_obj, slot = "data")  # Use 'data' for normalized counts
  
  # Ensure it's a matrix
  if (!is.matrix(assay_data)) {
    assay_data <- as.matrix(assay_data)
  }
  
  # Convert to a dataframe for manipulation
  expression_df <- as.data.frame(t(assay_data[genes_present, ]))  # Transpose for easier manipulation
  
  # Add cell type metadata
  expression_df$cell_type <- seurat_obj@meta.data$cellstate
  
  # Calculate average expression using dplyr summarise_at (mean)
  avg_expression <- expression_df %>%
    group_by(cell_type) %>%
    summarise_at(vars(all_of(genes_present)), mean, na.rm = TRUE) %>%
    pivot_longer(-cell_type, names_to = "gene", values_to = "avg_expression")
  
  # Z-score normalization across cell types for each gene
  avg_expression_long <- avg_expression %>%
    group_by(gene) %>%
    mutate(z_score = scale(avg_expression)) %>%
    ungroup()
  
  # Calculate the percentage of cells expressing each gene
  binary_expr <- assay_data[genes_present, ] > 0  # Convert to binary expression (expressed or not)
  binary_expr <- as.data.frame(t(binary_expr))  # Transpose for cell-wise rows
  
  # Add cell type information
  binary_expr$cell_type <- seurat_obj@meta.data$cellstate
  
  # Compute percentage of expressing cells per gene and cell type
  percentage_expr <- binary_expr %>%
    group_by(cell_type) %>%
    summarise(across(all_of(genes_present), ~ mean(.x) * 100, .names = "pct_{col}")) %>%
    pivot_longer(-cell_type, names_to = "gene", values_to = "percentage_expr") %>%
    mutate(gene = gsub("pct_", "", gene))  # Clean up gene names
  
  # Merge average expression and percentage expression data
  plot_data <- avg_expression_long %>%
    left_join(percentage_expr, by = c("gene", "cell_type")) %>%
    mutate(cell_type_dataset = paste(cell_type, dataset_name, sep = "_"))  # Combine cell type and dataset
  
  return(plot_data)
}

# Process each dataset manually
plot_data_Neftel <- process_dataset(obj_Neftel_Mal, "Neftel")
plot_data_Ebert <- process_dataset(obj_Ebert_Mal, "Ebert")
plot_data_LeBlanc <- process_dataset(obj_LeBlanc_Mal, "LeBlanc")
plot_data_Adelfattah <- process_dataset(obj_Abdelfattah_Mal, "Adelfattah")
plot_data_Chen <- process_dataset(obj_Chen_Mal, "Chen")
plot_data_Couturier <- process_dataset(obj_Couturier_Mal, "Couturier")

# Combine all datasets into a single dataframe for plotting
combined_plot_data <- bind_rows(
  plot_data_Couturier,
  plot_data_Chen,
  plot_data_Adelfattah,
  plot_data_LeBlanc,
  plot_data_Ebert,
  plot_data_Neftel
)

# Specify the order of genes for plotting
custom_gene_order <- c("DHODH", "PRMT5", "WDR77")

# Specify the order of cell annotations (cellstate) for plotting
custom_cellstate_order <- c("MES2", "MES1", "AC", "OPC", "NPC1", "NPC2")

# Ensure the 'gene' column respects the specified order
combined_plot_data$gene <- factor(combined_plot_data$gene, levels = custom_gene_order)

# Ensure the 'cellstate' in cell_type_dataset respects the specified order
combined_plot_data$cell_type_dataset <- factor(
  combined_plot_data$cell_type_dataset,
  levels = unlist(
    lapply(custom_cellstate_order, function(cellstate) {
      grep(paste0("^", cellstate, "_"), unique(combined_plot_data$cell_type_dataset), value = TRUE)
    })
  )
)

# Plot the combined dot plot with specified orders
ggplot(combined_plot_data, aes(x = gene, y = cell_type_dataset)) +
  geom_point(aes(size = percentage_expr, color = z_score)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Gene", y = "Cell Type (Dataset)", size = "Percentage Expressing", color = "Z-score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## permutation analysis to identify significantly enriched genes in a cellstate\

# Required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)

process_seurat_obj <- function(seurat_obj, Genesofinterest, num_permutations = 10000) {
  # Step 1: Extract the expression matrix and metadata
  expr_matrix <- GetAssayData(seurat_obj, slot = "data")
  metadata <- seurat_obj@meta.data
  
  # Step 2: Check if genes of interest are in the expression matrix
  valid_genes <- Genesofinterest %in% rownames(expr_matrix)
  if (!all(valid_genes)) {
    warning("Some genes in Genesofinterest are not found in the expression matrix. Proceeding with valid genes.")
    Genesofinterest <- Genesofinterest[valid_genes]
  }
  
  # Step 3: Filter the expression matrix for the genes of interest
  expr_matrix <- expr_matrix[Genesofinterest, , drop = FALSE]
  
  # Step 4: Get unique cell states
  cellstates <- unique(metadata$cellstate)
  
  # Step 5: Initialize a list to store average expressions
  avg_expr_list <- list()
  
  # Step 6: Iterate through each cell state
  for (cellstate in cellstates) {
    # Subset cells of the current cell state
    cell_indices <- which(metadata$cellstate == cellstate)
    cellstate_expr <- expr_matrix[, cell_indices]
    
    # Calculate the average expression for each gene in the cell state
    avg_expr_cellstate <- rowMeans(cellstate_expr)
    
    # Store the result in the list
    avg_expr_list[[cellstate]] <- avg_expr_cellstate
  }
  
  # Step 7: Convert the list to a dataframe
  average_expression <- do.call(rbind, avg_expr_list) %>%
    as.data.frame() %>%
    rownames_to_column(var = "cellstate") %>%
    gather(key = "gene", value = "avg_expression", -cellstate)
  
  # Step 8: Initialize the permutation results list
  permutation_results <- list()
  
  # Step 9: Perform 10,000 permutations
  set.seed(123) # For reproducibility
  
  for (i in 1:num_permutations) {
    # Randomly sample 1000 cells
    sampled_cells <- sample(colnames(seurat_obj), 1000, replace = TRUE)
    
    # Calculate average expression for the sampled cells
    permuted_expression <- expr_matrix[, sampled_cells] %>%
      as.data.frame() %>%
      rowMeans() %>%
      enframe(name = "gene", value = "perm_avg_expression")
    
    permutation_results[[i]] <- permuted_expression
  }
  
  # Step 10: Combine permutation results into a single dataframe
  permutation_df <- bind_rows(permutation_results, .id = "permutation") %>%
    group_by(gene) %>%
    summarize(perm_avg_expression = list(perm_avg_expression)) %>%
    ungroup()
  
  # Step 11: Calculate the proportion of times the random sample average is higher than each cell state average
  final_results <- average_expression %>%
    rowwise() %>%
    mutate(proportion_higher = mean(unlist(permutation_df$perm_avg_expression[permutation_df$gene == gene]) > avg_expression))
  
  # Step 12: Convert to dataframe and display the results
  final_results_df <- final_results %>% as.data.frame()
  
  # Step 13: Transform final_results_df to have genes as rows and cell states as columns
  final_results_wide <- final_results_df %>%
    select(gene, cellstate, proportion_higher) %>%
    spread(key = cellstate, value = proportion_higher)
  
  return(final_results_wide)
}

# Run the function with a specific Seurat dataset
final_results_wide <- process_seurat_obj(obj_Neftel_Mal, genes_of_interest)


############ Average expression of de novo pyrimidine synth pathway signature across Neftel cell states

# De novo pyrimidine synthesis gene signature
Genesofinterest <- c("CAD", "DHODH", "UMPS")

seurat_obj <- obj_Neftel_Mal #substitute with dataset of interest

# Load necessary library
library(dplyr)

# Define the z-score normalization function
z_score_normalize <- function(x) {
  Mean <- mean(x, na.rm = TRUE)
  SD <- sd(x, na.rm = TRUE)
  (x - Mean) / SD
}

# Extract the expression data for the Genesofinterest genes
Genesofinterest_features <- FetchData(object = seurat_obj, vars = Genesofinterest, slot = "data")

# Calculate the average expression of the Genesofinterest gene signature for each cell
Genesofinterest_avg_per_cell <- rowMeans(Genesofinterest_features, na.rm = TRUE)

# Add this average expression as metadata to the Seurat object
seurat_obj[["Genesofinterest_avg_Score"]] <- Genesofinterest_avg_per_cell

# Calculate the average Genesofinterest score per cell type using base R
Genesofinterest_avg_Score <- aggregate(Genesofinterest_avg_Score ~ cellstate, data = seurat_obj@meta.data, FUN = mean)

# Remove the "None" state assignment
Genesofinterest_avg_Score <- Genesofinterest_avg_Score[Genesofinterest_avg_Score$cellstate != "None", ]

# Apply the z-score normalization to the calculated average scores - do not do this step for permutation test below
Genesofinterest_avg_Score$Genesofinterest_avg_Score <- z_score_normalize(Genesofinterest_avg_Score$Genesofinterest_avg_Score)

# View the result
print(Genesofinterest_avg_Score)

## Permutation test for de novo pyrimidine synth pathway signature

# Required libraries
library(dplyr)
library(tidyr)
library(tibble)

# calculate the Genesofinterest_avg_Score (non z score normalized) for your Seurat object in your code block above

# Define the metadata column name
metadata_column <- "cellstate"

# Extract the metadata
metadata <- seurat_obj@meta.data

# Update the Seurat object to only include the filtered cells
seurat_obj <- subset(seurat_obj, cells = rownames(metadata))

# Get unique cell states after filtering
celltypes <- unique(metadata[[metadata_column]])

# Initialize a list to store average Genesofinterest scores
avg_Genesofinterest_list <- list()

# Iterate through each cell state
for (celltype in celltypes) {
  # Subset cells of the current cell state
  cell_indices <- which(metadata_filtered[[metadata_column]] == celltype)
  
  # Check if there are any cells of the current cell state
  if (length(cell_indices) > 0) {
    # Calculate the average Genesofinterest score for the current cell state
    avg_Genesofinterest_celltype <- mean(metadata_filtered$Genesofinterest_avg_Score[cell_indices], na.rm = TRUE)
    
    # Store the result in the list
    avg_Genesofinterest_list[[celltype]] <- avg_Genesofinterest_celltype
  } else {
    warning(paste("No cells found for cell state:", celltype))
  }
}

# Convert the list to a dataframe
average_Genesofinterest_scores <- data.frame(
  celltype = names(avg_Genesofinterest_list),
  avg_Genesofinterest_score = unlist(avg_Genesofinterest_list)
)

# Initialize the permutation results list
permutation_results <- list()

# Perform 10,000 permutations
set.seed(123) # For reproducibility
num_permutations <- 10000
for (i in 1:num_permutations) {
  # Randomly shuffle the Genesofinterest_avg_Score within the filtered metadata
  permuted_scores <- sample(metadata_filtered$Genesofinterest_avg_Score, replace = TRUE)
  
  # Calculate average Genesofinterest score for each cell type with permuted data
  perm_avg_Genesofinterest_list <- list()
  for (celltype in celltypes) {
    # Subset cells of the current cell state
    cell_indices <- which(metadata_filtered[[metadata_column]] == celltype)
    
    # Calculate the average permuted Genesofinterest score
    if (length(cell_indices) > 0) {
      perm_avg_Genesofinterest_list[[celltype]] <- mean(permuted_scores[cell_indices], na.rm = TRUE)
    }
  }
  
  # Store the permutation results
  permutation_results[[i]] <- perm_avg_Genesofinterest_list
}

# Convert permutation results into a dataframe
permutation_df <- bind_rows(lapply(permutation_results, function(x) as.data.frame(t(unlist(x)))), .id = "permutation")

# Calculate the proportion of times the random sample average is higher than each cell state average
final_results <- average_Genesofinterest_scores %>%
  rowwise() %>%
  mutate(proportion_higher = mean(permutation_df[[celltype]] > avg_Genesofinterest_score, na.rm = TRUE))

# Convert to dataframe and display the results
final_results_df <- final_results %>% as.data.frame()

# Transform final_results_df to have cell states as rows and the proportion_higher column
final_results_wide <- final_results_df %>%
  select(celltype, proportion_higher)

# Print the final_results_wide to check the results
print(final_results_wide)
