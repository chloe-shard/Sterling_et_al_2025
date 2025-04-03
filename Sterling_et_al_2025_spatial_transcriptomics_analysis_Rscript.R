#################### Spatial transcriptomics analysis #####################
###########################################################################

##load spatial files

# Required library
library(Seurat)

# Define the path where the Seurat objects are saved
path_to_files <- "C:/Users/shardcl/OneDrive - University of South Australia/NHMRC Ideas 2021 GNT 2013180/Spatial Transcriptomics/RCTD_plots"

# List of Seurat object names
seurat_names <- c("T243_RCTD_all_subpop_cnv", "T241_RCTD_all_subpop_cnv", "T242_RCTD_all_subpop_cnv", "T248_RCTD_all_subpop_cnv", 
                  "T251_RCTD_all_subpop_cnv", "T255_RCTD_all_subpop_cnv", "T256_RCTD_all_subpop_cnv", "T259_RCTD_all_subpop_cnv", 
                  "T260_RCTD_all_subpop_cnv", "T262_RCTD_all_subpop_cnv", "T265_RCTD_all_subpop_cnv", "T266_RCTD_all_subpop_cnv", 
                  "T268_RCTD_all_subpop_cnv", "T269_RCTD_all_subpop_cnv", "T270_RCTD_all_subpop_cnv", "T275_RCTD_all_subpop_cnv", 
                  "T296_RCTD_all_subpop_cnv", "T304_RCTD_all_subpop_cnv", "T313_RCTD_all_subpop_cnv", "T334_RCTD_all_subpop_cnv", 
                  "GBM_MGH258_RCTD_all_subpop_cnv", "GBM_ZH881inf_RCTD_all_subpop_cnv", "GBM_ZH881T1_RCTD_all_subpop_cnv", 
                  "GBM_ZH916bulk_RCTD_all_subpop_cnv", "GBM_ZH916inf_RCTD_all_subpop_cnv", "GBM_ZH916T1_RCTD_all_subpop_cnv", 
                  "GBM_ZH1007inf_RCTD_all_subpop_cnv", "GBM_ZH1007nec_RCTD_all_subpop_cnv", "GBM_ZH1019inf_RCTD_all_subpop_cnv", 
                  "GBM_ZH1019T1_RCTD_all_subpop_cnv", "GBM_ZH8811Abulk_RCTD_all_subpop_cnv", "GBM_ZH8811Bbulk_RCTD_all_subpop_cnv", 
                  "GBM_ZH8812bulk_RCTD_all_subpop_cnv")

# Load each Seurat object if the file exists
for (seurat_name in seurat_names) {
  file_path <- file.path(path_to_files, paste0(seurat_name, ".rds"))
  if (file.exists(file_path)) {
    assign(seurat_name, readRDS(file_path))
    message(paste("Loaded:", seurat_name))
  } else {
    warning(paste("File not found:", file_path))
  }
}

##### SCALOP annotation of tumour regions

# Required libraries
library(readr)
library(scalop)

seurat_names <- c("T243_RCTD_all_subpop_cnv", "T241_RCTD_all_subpop_cnv", "T242_RCTD_all_subpop_cnv", "T248_RCTD_all_subpop_cnv", 
                  "T251_RCTD_all_subpop_cnv", "T255_RCTD_all_subpop_cnv", "T256_RCTD_all_subpop_cnv", "T259_RCTD_all_subpop_cnv", 
                  "T260_RCTD_all_subpop_cnv", "T262_RCTD_all_subpop_cnv", "T265_RCTD_all_subpop_cnv", "T266_RCTD_all_subpop_cnv", 
                  "T268_RCTD_all_subpop_cnv", "T269_RCTD_all_subpop_cnv", "T270_RCTD_all_subpop_cnv", "T275_RCTD_all_subpop_cnv", 
                  "T296_RCTD_all_subpop_cnv", "T304_RCTD_all_subpop_cnv", "T313_RCTD_all_subpop_cnv", "T334_RCTD_all_subpop_cnv", 
                  "GBM_MGH258_RCTD_all_subpop_cnv", "GBM_ZH881inf_RCTD_all_subpop_cnv", "GBM_ZH881T1_RCTD_all_subpop_cnv", 
                  "GBM_ZH916bulk_RCTD_all_subpop_cnv", "GBM_ZH916inf_RCTD_all_subpop_cnv", "GBM_ZH916T1_RCTD_all_subpop_cnv", 
                  "GBM_ZH1007inf_RCTD_all_subpop_cnv", "GBM_ZH1007nec_RCTD_all_subpop_cnv", "GBM_ZH1019inf_RCTD_all_subpop_cnv", 
                  "GBM_ZH1019T1_RCTD_all_subpop_cnv", "GBM_ZH8811Abulk_RCTD_all_subpop_cnv", "GBM_ZH8811Bbulk_RCTD_all_subpop_cnv", 
                  "GBM_ZH8812bulk_RCTD_all_subpop_cnv")

#load region signatures
Region_signatures <- read_csv("file location of Region_signatures.csv")

LE_genes <- Region_signatures[["LE"]]
CTmvp_genes <- Region_signatures[["CTmvp"]]
CTpan_genes <- Region_signatures[["CTpan"]]
CT_genes <- Region_signatures[["CT"]]

# List of signatures
sigs <- list(CT_genes, CTmvp_genes, CTpan_genes, LE_genes)
names(sigs) <- c("CT_genes", "CTmvp_genes", "CTpan_genes", "LE_genes")

# Loop through each sample in seurat_names
for (sample_name in seurat_names) {
  # Load the Seurat object if not already in the environment
  seurat_obj <- get(sample_name)
  
  # Basic Scoring of Matrix by Gene sigs
  m <- as.matrix(GetAssayData(object = seurat_obj[["Spatial"]], slot = "data"))
  
  m_scores <- baseScores(m, sigs, conserved.genes = 0.7)
  
  # Score a Matrix by Gene sigs (Signatures)
  m_sigscores <- sigScores(m, sigs, groups = NULL, center.rows = TRUE,
                           center = TRUE,
                           expr.center = TRUE,
                           expr.bin.m = NULL,
                           expr.bins = NULL,
                           expr.sigs = NULL,
                           expr.nbin = 30,
                           expr.binsize = 100,
                           conserved.genes = 0.7,
                           replace = FALSE
  )
  
  vector_max_sig <- maxcol_strict(m_sigscores, min = NULL, diff = NULL, splitByCol = FALSE)
  
  highest_Sig_Column <- as.data.frame(colnames(m_sigscores)[max.col(m_sigscores, ties.method = "first")])
  row.names(highest_Sig_Column) <- row.names(m_sigscores)
  highest_Sig_Column$Spot <- rownames(highest_Sig_Column)
  
  # Match the order of Seurat barcodes with your data
  myBarcode <- rownames(seurat_obj@meta.data)
  highest_Sig_Column_2 <- highest_Sig_Column[match(myBarcode, highest_Sig_Column$Spot), ]
  
  # Add the calculated values to the Seurat object's metadata
  seurat_obj$Tumour_region_SpotIDs <- highest_Sig_Column_2[, 1] # Adjust the column index if needed
  
  # Optionally, save or assign the modified Seurat object back to the environment
  assign(sample_name, seurat_obj)
  
  # Progress update (optional)
  print(paste("Processed:", sample_name))
}

# Examine region annotation
table(T266_RCTD_all_subpop_cnv$Tumour_region_SpotIDs)

############# Calculate aggregation score across tissues for tumour regions #############

# Required library
library(Seurat)

# List of Seurat objects
seurat_objects <- list(
  T243_RCTD_all_subpop_cnv = T243_RCTD_all_subpop_cnv, 
  T241_RCTD_all_subpop_cnv = T241_RCTD_all_subpop_cnv, 
  T242_RCTD_all_subpop_cnv = T242_RCTD_all_subpop_cnv, 
  T248_RCTD_all_subpop_cnv = T248_RCTD_all_subpop_cnv, 
  T251_RCTD_all_subpop_cnv = T251_RCTD_all_subpop_cnv, 
  T255_RCTD_all_subpop_cnv = T255_RCTD_all_subpop_cnv, 
  T256_RCTD_all_subpop_cnv = T256_RCTD_all_subpop_cnv, 
  T259_RCTD_all_subpop_cnv = T259_RCTD_all_subpop_cnv, 
  T260_RCTD_all_subpop_cnv = T260_RCTD_all_subpop_cnv, 
  T262_RCTD_all_subpop_cnv = T262_RCTD_all_subpop_cnv, 
  T265_RCTD_all_subpop_cnv = T265_RCTD_all_subpop_cnv, 
  T266_RCTD_all_subpop_cnv = T266_RCTD_all_subpop_cnv, 
  T268_RCTD_all_subpop_cnv = T268_RCTD_all_subpop_cnv, 
  T269_RCTD_all_subpop_cnv = T269_RCTD_all_subpop_cnv, 
  T270_RCTD_all_subpop_cnv = T270_RCTD_all_subpop_cnv, 
  T275_RCTD_all_subpop_cnv = T275_RCTD_all_subpop_cnv, 
  T296_RCTD_all_subpop_cnv = T296_RCTD_all_subpop_cnv, 
  T304_RCTD_all_subpop_cnv = T304_RCTD_all_subpop_cnv, 
  T313_RCTD_all_subpop_cnv = T313_RCTD_all_subpop_cnv, 
  T334_RCTD_all_subpop_cnv = T334_RCTD_all_subpop_cnv, 
  GBM_MGH258_RCTD_all_subpop_cnv = GBM_MGH258_RCTD_all_subpop_cnv, 
  GBM_ZH881inf_RCTD_all_subpop_cnv = GBM_ZH881inf_RCTD_all_subpop_cnv, 
  GBM_ZH881T1_RCTD_all_subpop_cnv = GBM_ZH881T1_RCTD_all_subpop_cnv, 
  GBM_ZH916bulk_RCTD_all_subpop_cnv = GBM_ZH916bulk_RCTD_all_subpop_cnv, 
  GBM_ZH916inf_RCTD_all_subpop_cnv = GBM_ZH916inf_RCTD_all_subpop_cnv, 
  GBM_ZH916T1_RCTD_all_subpop_cnv = GBM_ZH916T1_RCTD_all_subpop_cnv, 
  GBM_ZH1007inf_RCTD_all_subpop_cnv = GBM_ZH1007inf_RCTD_all_subpop_cnv, 
  GBM_ZH1007nec_RCTD_all_subpop_cnv = GBM_ZH1007nec_RCTD_all_subpop_cnv, 
  GBM_ZH1019inf_RCTD_all_subpop_cnv = GBM_ZH1019inf_RCTD_all_subpop_cnv, 
  GBM_ZH1019T1_RCTD_all_subpop_cnv = GBM_ZH1019T1_RCTD_all_subpop_cnv, 
  GBM_ZH8811Abulk_RCTD_all_subpop_cnv = GBM_ZH8811Abulk_RCTD_all_subpop_cnv, 
  GBM_ZH8811Bbulk_RCTD_all_subpop_cnv = GBM_ZH8811Bbulk_RCTD_all_subpop_cnv, 
  GBM_ZH8812bulk_RCTD_all_subpop_cnv = GBM_ZH8812bulk_RCTD_all_subpop_cnv
)


# Specify the metadata column and spot annotation to use
metadata_column <- "Tumour_region_SpotIDs"
spot_annotation <- "LE_genes"  # Specify the annotation you want to process

# Initialize a dataframe to store the results
results <- data.frame(SeuratObject = character(), AverageFractionPositive = numeric(), stringsAsFactors = FALSE)

# Function to process each Seurat object
process_spatial_object <- function(spatial_obj, spot_annotation) {
  # Check if the specified metadata column is present
  if (!metadata_column %in% colnames(spatial_obj@meta.data)) {
    warning(paste("Annotation", metadata_column, "not found in object:", deparse(substitute(spatial_obj))))
    return(NA)
  }
  
  # Get the tissue coordinates
  coords <- GetTissueCoordinates(spatial_obj)
  
  # Identify positive spots for the specific spot annotation
  positive_spots <- WhichCells(spatial_obj, cells = rownames(spatial_obj@meta.data[spatial_obj@meta.data[[metadata_column]] == spot_annotation,]))
  
  # Check if there are any positive spots
  if (length(positive_spots) == 0) {
    warning(paste("No positive spots found for annotation", spot_annotation, "in object:", deparse(substitute(spatial_obj))))
    return(NA)
  }
  
  # Get coordinates of the positive spots
  positive_coords <- coords[positive_spots, ]
  
  # Create an empty list to store distance matrices
  distance_matrices <- list()
  
  # Function to calculate distances
  calculate_distances <- function(spot_coords, all_coords) {
    # Ensure spot_coords is a vector
    spot_vector <- as.numeric(spot_coords)
    # Calculate Euclidean distances
    distances <- sqrt(rowSums((t(t(all_coords) - spot_vector)) ^ 2))
    return(distances)
  }
  
  # Loop over each positive spot and calculate distances to all other spots
  for (spot in positive_spots) {
    # Get the coordinates of the positive spot
    spot_coords <- coords[spot, , drop = FALSE]
    
    # Calculate distances
    distances <- calculate_distances(spot_coords, coords)
    
    # Create a data frame with SpotID and distances, ordered by distance
    distance_df <- data.frame(SpotID = rownames(coords), Distance = distances)
    distance_df <- distance_df[order(distance_df$Distance), ]
    
    # Store the distance data frame in the list
    distance_matrices[[spot]] <- distance_df
  }
  
  # Calculate the distance threshold based on the value halfway between the average distance of the 6th and 7th closest neighbors
  sixth_distances <- sapply(distance_matrices, function(df) df$Distance[7])
  seventh_distances <- sapply(distance_matrices, function(df) df$Distance[8])
  distance_threshold <- mean(sixth_distances + seventh_distances) / 2
  
  # Create an empty matrix with rows as positive spots and columns as 6 neighbors
  neighbor_matrix <- matrix(NA, nrow = length(positive_spots), ncol = 6)
  rownames(neighbor_matrix) <- positive_spots
  colnames(neighbor_matrix) <- paste0("Neighbor_", 1:6)
  
  # Create a vector to store all neighbor SpotIDs
  all_neighbors <- c()
  
  # Loop over each positive spot to fill the neighbor matrix
  for (i in seq_along(positive_spots)) {
    spot <- positive_spots[i]
    distance_df <- distance_matrices[[spot]]
    
    # Filter neighbors within the distance threshold (excluding the spot itself)
    neighbors <- distance_df$SpotID[distance_df$Distance <= distance_threshold & distance_df$SpotID != spot]
    
    # Fill the matrix with the first 6 neighbors
    neighbor_matrix[i, 1:min(length(neighbors), 6)] <- neighbors[1:min(length(neighbors), 6)]
    
    # Collect all neighbor SpotIDs
    all_neighbors <- unique(c(all_neighbors, neighbors))
  }
  
  # Replace SpotIDs with the annotations in the specified metadata column
  neighbor_matrix_annotated <- neighbor_matrix
  for (i in 1:nrow(neighbor_matrix)) {
    for (j in 1:ncol(neighbor_matrix)) {
      spot_id <- neighbor_matrix[i, j]
      if (!is.na(spot_id)) {
        neighbor_matrix_annotated[i, j] <- ifelse(spatial_obj@meta.data[spot_id, metadata_column] == spot_annotation, "positive", "negative")
      }
    }
  }
  
  # Calculate the fraction of "positive" annotations for each row
  positive_fractions <- apply(neighbor_matrix_annotated, 1, function(row) {
    positive_count <- sum(row == "positive", na.rm = TRUE)
    total_count <- sum(!is.na(row))
    fraction_positive <- positive_count / total_count
    return(fraction_positive)
  })
  
  # Calculate the average fraction score across all rows
  average_fraction_positive <- mean(positive_fractions, na.rm = TRUE)
  
  return(average_fraction_positive)
}

# Process each Seurat object for the specified spot annotation and store the results
for (obj_name in names(seurat_objects)) {
  spatial_obj <- seurat_objects[[obj_name]]
  if (metadata_column %in% colnames(spatial_obj@meta.data)) {
    avg_fraction_positive <- process_spatial_object(spatial_obj, spot_annotation)
    if (!is.na(avg_fraction_positive)) {
      results <- rbind(results, data.frame(SeuratObject = obj_name, AverageFractionPositive = avg_fraction_positive))
    }
  } else {
    warning(paste("Annotation", metadata_column, "not found in object:", obj_name))
  }
}

# View the results
print(results)

############### Average expression of genes in tumour regions (within vs outside)
#note: region enriched tissues are tiisues with an aggregation score >0.4 (calulated in code bloack above)

#LE enriched
seurat_objects <- list(
  T243_RCTD_all_subpop_cnv = T243_RCTD_all_subpop_cnv,
  T241_RCTD_all_subpop_cnv = T241_RCTD_all_subpop_cnv,
  T242_RCTD_all_subpop_cnv = T242_RCTD_all_subpop_cnv,
  T256_RCTD_all_subpop_cnv = T256_RCTD_all_subpop_cnv,
  T260_RCTD_all_subpop_cnv = T260_RCTD_all_subpop_cnv,
  T262_RCTD_all_subpop_cnv = T262_RCTD_all_subpop_cnv,
  T266_RCTD_all_subpop_cnv = T266_RCTD_all_subpop_cnv,
  T269_RCTD_all_subpop_cnv = T269_RCTD_all_subpop_cnv,
  T270_RCTD_all_subpop_cnv = T270_RCTD_all_subpop_cnv,
  T275_RCTD_all_subpop_cnv = T275_RCTD_all_subpop_cnv,
  T296_RCTD_all_subpop_cnv = T296_RCTD_all_subpop_cnv,
  T304_RCTD_all_subpop_cnv = T304_RCTD_all_subpop_cnv,
  T334_RCTD_all_subpop_cnv = T334_RCTD_all_subpop_cnv,
  GBM_ZH881inf_RCTD_all_subpop_cnv = GBM_ZH881inf_RCTD_all_subpop_cnv,
  GBM_ZH881T1_RCTD_all_subpop_cnv = GBM_ZH881T1_RCTD_all_subpop_cnv,
  GBM_ZH916bulk_RCTD_all_subpop_cnv = GBM_ZH916bulk_RCTD_all_subpop_cnv,
  GBM_ZH916inf_RCTD_all_subpop_cnv = GBM_ZH916inf_RCTD_all_subpop_cnv,
  GBM_ZH1019inf_RCTD_all_subpop_cnv = GBM_ZH1019inf_RCTD_all_subpop_cnv,
  GBM_ZH8811Abulk_RCTD_all_subpop_cnv = GBM_ZH8811Abulk_RCTD_all_subpop_cnv,
  GBM_ZH8812bulk_RCTD_all_subpop_cnv = GBM_ZH8812bulk_RCTD_all_subpop_cnv
)

#CT enriched
seurat_objects <- list(
  T243_RCTD_all_subpop_cnv = T243_RCTD_all_subpop_cnv,
  T242_RCTD_all_subpop_cnv = T242_RCTD_all_subpop_cnv,
  T248_RCTD_all_subpop_cnv = T248_RCTD_all_subpop_cnv,
  T251_RCTD_all_subpop_cnv = T251_RCTD_all_subpop_cnv,
  T255_RCTD_all_subpop_cnv = T255_RCTD_all_subpop_cnv,
  T256_RCTD_all_subpop_cnv = T256_RCTD_all_subpop_cnv,
  T259_RCTD_all_subpop_cnv = T259_RCTD_all_subpop_cnv,
  T260_RCTD_all_subpop_cnv = T260_RCTD_all_subpop_cnv,
  T262_RCTD_all_subpop_cnv = T262_RCTD_all_subpop_cnv,
  T266_RCTD_all_subpop_cnv = T266_RCTD_all_subpop_cnv,
  T268_RCTD_all_subpop_cnv = T268_RCTD_all_subpop_cnv,
  T269_RCTD_all_subpop_cnv = T269_RCTD_all_subpop_cnv,
  T275_RCTD_all_subpop_cnv = T275_RCTD_all_subpop_cnv,
  T304_RCTD_all_subpop_cnv = T304_RCTD_all_subpop_cnv,
  T313_RCTD_all_subpop_cnv = T313_RCTD_all_subpop_cnv,
  T334_RCTD_all_subpop_cnv = T334_RCTD_all_subpop_cnv,
  GBM_MGH258_RCTD_all_subpop_cnv = GBM_MGH258_RCTD_all_subpop_cnv,
  GBM_ZH881inf_RCTD_all_subpop_cnv = GBM_ZH881inf_RCTD_all_subpop_cnv,
  GBM_ZH881T1_RCTD_all_subpop_cnv = GBM_ZH881T1_RCTD_all_subpop_cnv,
  GBM_ZH916bulk_RCTD_all_subpop_cnv = GBM_ZH916bulk_RCTD_all_subpop_cnv,
  GBM_ZH916inf_RCTD_all_subpop_cnv = GBM_ZH916inf_RCTD_all_subpop_cnv,
  GBM_ZH916T1_RCTD_all_subpop_cnv = GBM_ZH916T1_RCTD_all_subpop_cnv,
  GBM_ZH1007inf_RCTD_all_subpop_cnv = GBM_ZH1007inf_RCTD_all_subpop_cnv,
  GBM_ZH1007nec_RCTD_all_subpop_cnv = GBM_ZH1007nec_RCTD_all_subpop_cnv,
  GBM_ZH1019inf_RCTD_all_subpop_cnv = GBM_ZH1019inf_RCTD_all_subpop_cnv,
  GBM_ZH1019T1_RCTD_all_subpop_cnv = GBM_ZH1019T1_RCTD_all_subpop_cnv,
  GBM_ZH8811Abulk_RCTD_all_subpop_cnv = GBM_ZH8811Abulk_RCTD_all_subpop_cnv,
  GBM_ZH8811Bbulk_RCTD_all_subpop_cnv = GBM_ZH8811Bbulk_RCTD_all_subpop_cnv,
  GBM_ZH8812bulk_RCTD_all_subpop_cnv = GBM_ZH8812bulk_RCTD_all_subpop_cnv
)

#CTmvp enriched
seurat_objects <- list(
  T248_RCTD_all_subpop_cnv = T248_RCTD_all_subpop_cnv,
  T313_RCTD_all_subpop_cnv = T313_RCTD_all_subpop_cnv,
  GBM_ZH881T1_RCTD_all_subpop_cnv = GBM_ZH881T1_RCTD_all_subpop_cnv,
  GBM_ZH916bulk_RCTD_all_subpop_cnv = GBM_ZH916bulk_RCTD_all_subpop_cnv,
  GBM_ZH916T1_RCTD_all_subpop_cnv = GBM_ZH916T1_RCTD_all_subpop_cnv,
  GBM_ZH1007inf_RCTD_all_subpop_cnv = GBM_ZH1007inf_RCTD_all_subpop_cnv,
  GBM_ZH1007nec_RCTD_all_subpop_cnv = GBM_ZH1007nec_RCTD_all_subpop_cnv
)

#CTpan enriched
seurat_objects <- list(
  T243_RCTD_all_subpop_cnv = T243_RCTD_all_subpop_cnv,
  T248_RCTD_all_subpop_cnv = T248_RCTD_all_subpop_cnv,
  T255_RCTD_all_subpop_cnv = T255_RCTD_all_subpop_cnv,
  T259_RCTD_all_subpop_cnv = T259_RCTD_all_subpop_cnv,
  T260_RCTD_all_subpop_cnv = T260_RCTD_all_subpop_cnv,
  T265_RCTD_all_subpop_cnv = T265_RCTD_all_subpop_cnv,
  T269_RCTD_all_subpop_cnv = T269_RCTD_all_subpop_cnv,
  T275_RCTD_all_subpop_cnv = T275_RCTD_all_subpop_cnv,
  T304_RCTD_all_subpop_cnv = T304_RCTD_all_subpop_cnv,
  T313_RCTD_all_subpop_cnv = T313_RCTD_all_subpop_cnv,
  T334_RCTD_all_subpop_cnv = T334_RCTD_all_subpop_cnv,
  GBM_ZH881T1_RCTD_all_subpop_cnv = GBM_ZH881T1_RCTD_all_subpop_cnv,
  GBM_ZH916bulk_RCTD_all_subpop_cnv = GBM_ZH916bulk_RCTD_all_subpop_cnv,
  GBM_ZH916T1_RCTD_all_subpop_cnv = GBM_ZH916T1_RCTD_all_subpop_cnv,
  GBM_ZH1007inf_RCTD_all_subpop_cnv = GBM_ZH1007inf_RCTD_all_subpop_cnv,
  GBM_ZH1007nec_RCTD_all_subpop_cnv = GBM_ZH1007nec_RCTD_all_subpop_cnv,
  GBM_ZH1019inf_RCTD_all_subpop_cnv = GBM_ZH1019inf_RCTD_all_subpop_cnv,
  GBM_ZH1019T1_RCTD_all_subpop_cnv = GBM_ZH1019T1_RCTD_all_subpop_cnv,
  GBM_ZH8811Abulk_RCTD_all_subpop_cnv = GBM_ZH8811Abulk_RCTD_all_subpop_cnv,
  GBM_ZH8811Bbulk_RCTD_all_subpop_cnv = GBM_ZH8811Bbulk_RCTD_all_subpop_cnv
)


# Define the gene of interest
gene_of_interest <- "DHODH"  # Replace with the actual gene name

# Define the region of interest in Tumour_region_SpotIDs
region_of_interest <- "CTmvp_genes"  #(make sure to select corresponding Seurat object list above)

# Initialize an empty dataframe to store results
results <- data.frame(
  Sample = character(),
  Region = character(),
  AvgExpression_Region = numeric(),
  AvgExpression_non_Region = numeric(),
  log2FC = numeric(),
  stringsAsFactors = FALSE
)

# Set a seed for reproducibility
set.seed(1234)

# Iterate over each Seurat object
for (sample_name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[sample_name]]
  
  # Ensure the gene exists in the dataset
  if (!(gene_of_interest %in% rownames(seurat_obj))) {
    message(paste("Skipping", sample_name, "- Gene not found"))
    next
  }
  
  # Get all spots labeled as the chosen region of interest
  region_spots <- WhichCells(seurat_obj, expression = Tumour_region_SpotIDs == region_of_interest)
  
  # Get all non-selected spots
  non_region_spots <- setdiff(Cells(seurat_obj), region_spots)
  
  # Check if there are any region-specific spots
  if (length(region_spots) == 0 || length(non_region_spots) == 0) {
    message(paste("Skipping", sample_name, "- No spots found for", region_of_interest))
    next
  }
  
  # Determine the number of spots to use (minimum of region or non-region)
  num_spots <- min(length(region_spots), length(non_region_spots))
  
  # Set seed before sampling to ensure reproducibility
  set.seed(1234)
  selected_region_spots <- sample(region_spots, num_spots)
  
  set.seed(1234)
  selected_non_region_spots <- sample(non_region_spots, num_spots)
  
  # Extract expression data for the gene in selected spots
  gene_expression_region <- FetchData(seurat_obj, vars = gene_of_interest, cells = selected_region_spots)
  avg_expression_region <- mean(gene_expression_region[[gene_of_interest]], na.rm = TRUE)
  
  # Extract expression data for the gene in randomly selected non-region spots
  gene_expression_non_region <- FetchData(seurat_obj, vars = gene_of_interest, cells = selected_non_region_spots)
  avg_expression_non_region <- mean(gene_expression_non_region[[gene_of_interest]], na.rm = TRUE)
  
  # Calculate log2 Fold Change (log2FC)
  log2_FC <- log2((avg_expression_region + 1) / (avg_expression_non_region + 1))
  
  # Store results
  results <- rbind(
    results, 
    data.frame(
      Sample = sample_name, 
      Region = region_of_interest,
      AvgExpression_Region = avg_expression_region, 
      AvgExpression_non_Region = avg_expression_non_region,
      log2FC = log2_FC
    )
  )
}

# Display results
print(results)

################ De novo pyrimidine synthesis signature exp across regions

de_novo_genes <- c("CAD", "DHODH", "UMPS")

# List of Seurat objects
seurat_objects <- list(
  T243_RCTD_all_subpop_cnv = T243_RCTD_all_subpop_cnv, 
  T241_RCTD_all_subpop_cnv = T241_RCTD_all_subpop_cnv, 
  T242_RCTD_all_subpop_cnv = T242_RCTD_all_subpop_cnv, 
  T248_RCTD_all_subpop_cnv = T248_RCTD_all_subpop_cnv, 
  T251_RCTD_all_subpop_cnv = T251_RCTD_all_subpop_cnv, 
  T255_RCTD_all_subpop_cnv = T255_RCTD_all_subpop_cnv, 
  T256_RCTD_all_subpop_cnv = T256_RCTD_all_subpop_cnv, 
  T259_RCTD_all_subpop_cnv = T259_RCTD_all_subpop_cnv, 
  T260_RCTD_all_subpop_cnv = T260_RCTD_all_subpop_cnv, 
  T262_RCTD_all_subpop_cnv = T262_RCTD_all_subpop_cnv, 
  T265_RCTD_all_subpop_cnv = T265_RCTD_all_subpop_cnv, 
  T266_RCTD_all_subpop_cnv = T266_RCTD_all_subpop_cnv, 
  T268_RCTD_all_subpop_cnv = T268_RCTD_all_subpop_cnv, 
  T269_RCTD_all_subpop_cnv = T269_RCTD_all_subpop_cnv, 
  T270_RCTD_all_subpop_cnv = T270_RCTD_all_subpop_cnv, 
  T275_RCTD_all_subpop_cnv = T275_RCTD_all_subpop_cnv, 
  T296_RCTD_all_subpop_cnv = T296_RCTD_all_subpop_cnv, 
  T304_RCTD_all_subpop_cnv = T304_RCTD_all_subpop_cnv, 
  T313_RCTD_all_subpop_cnv = T313_RCTD_all_subpop_cnv, 
  T334_RCTD_all_subpop_cnv = T334_RCTD_all_subpop_cnv, 
  GBM_MGH258_RCTD_all_subpop_cnv = GBM_MGH258_RCTD_all_subpop_cnv, 
  GBM_ZH881inf_RCTD_all_subpop_cnv = GBM_ZH881inf_RCTD_all_subpop_cnv, 
  GBM_ZH881T1_RCTD_all_subpop_cnv = GBM_ZH881T1_RCTD_all_subpop_cnv, 
  GBM_ZH916bulk_RCTD_all_subpop_cnv = GBM_ZH916bulk_RCTD_all_subpop_cnv, 
  GBM_ZH916inf_RCTD_all_subpop_cnv = GBM_ZH916inf_RCTD_all_subpop_cnv, 
  GBM_ZH916T1_RCTD_all_subpop_cnv = GBM_ZH916T1_RCTD_all_subpop_cnv, 
  GBM_ZH1007inf_RCTD_all_subpop_cnv = GBM_ZH1007inf_RCTD_all_subpop_cnv, 
  GBM_ZH1007nec_RCTD_all_subpop_cnv = GBM_ZH1007nec_RCTD_all_subpop_cnv, 
  GBM_ZH1019inf_RCTD_all_subpop_cnv = GBM_ZH1019inf_RCTD_all_subpop_cnv, 
  GBM_ZH1019T1_RCTD_all_subpop_cnv = GBM_ZH1019T1_RCTD_all_subpop_cnv, 
  GBM_ZH8811Abulk_RCTD_all_subpop_cnv = GBM_ZH8811Abulk_RCTD_all_subpop_cnv, 
  GBM_ZH8811Bbulk_RCTD_all_subpop_cnv = GBM_ZH8811Bbulk_RCTD_all_subpop_cnv, 
  GBM_ZH8812bulk_RCTD_all_subpop_cnv = GBM_ZH8812bulk_RCTD_all_subpop_cnv
)


# Named list of gene signatures
gene_signatures <- list(
  de_novo = de_novo_genes,
  salvage = salvage_genes,
  PRMT5_WDR77 = PRMT5_WDR77_complex
)


# Function to process each Seurat object
process_seurat_object <- function(seurat_obj, gene_signatures) {
  
  # Create an empty list to store the average scores
  avg_scores_list <- list()
  
  for (signature_name in names(gene_signatures)) {
    
    # Fetch data for genes present in the current Seurat object
    gene_list <- gene_signatures[[signature_name]]
    genes_in_obj <- intersect(gene_list, rownames(seurat_obj))
    
    if (length(genes_in_obj) > 0) {
      # Fetch data, ensuring cells are rows and genes are columns
      genes_features <- FetchData(object = seurat_obj, vars = genes_in_obj, slot = "data")
      
      # Calculate the average score per cell (spot) across the selected genes
      avg_score <- rowMeans(genes_features, na.rm = TRUE)  # Correct calculation using rowMeans
      
      # Add the average scores to the list
      avg_scores_list[[signature_name]] <- avg_score
    } else {
      # Add a column of NAs if no genes are found in the Seurat object
      avg_scores_list[[signature_name]] <- rep(NA, ncol(seurat_obj))
    }
  }
  
  # Combine all the average scores into a single dataframe
  Tumour_region <- as.data.frame(avg_scores_list)
  
  # Ensure the row names of Tumour_region match the cell names in the Seurat object
  rownames(Tumour_region) <- colnames(seurat_obj)
  
  # Diagnostic print statements to check alignment
  message("Number of spots (cells) in Seurat object: ", ncol(seurat_obj))
  message("Number of rows in Tumour_region: ", nrow(Tumour_region))
  message("Column names in Tumour_region: ", paste(colnames(Tumour_region), collapse = ", "))
  message("Adding metadata to Seurat object...")
  
  # Add metadata to the Seurat object
  seurat_obj <- AddMetaData(seurat_obj, Tumour_region)
  
  # Verify if metadata was added
  if (all(colnames(Tumour_region) %in% names(seurat_obj@meta.data))) {
    message("Successfully added metadata columns: ", paste(colnames(Tumour_region), collapse = ", "))
  } else {
    message("Failed to add some or all metadata columns.")
  }
  
  return(seurat_obj)
}

# Process each Seurat object in the list and update the list
seurat_objects <- lapply(names(seurat_objects), function(obj_name) {
  # Output progress message
  message("Processing ", obj_name)
  
  # Process and update the Seurat object
  seurat_objects[[obj_name]] <- process_seurat_object(seurat_objects[[obj_name]], gene_signatures)
  
  # Reassign the modified object to the global environment
  assign(obj_name, seurat_objects[[obj_name]], envir = .GlobalEnv)
  
  return(seurat_objects[[obj_name]])
})

# Check if the metadata columns were added correctly
head(T243_RCTD_all_subpop_cnv@meta.data)

# Plot the Genesofinterest score (optional)
SpatialFeaturePlot(T260_RCTD_all_subpop_cnv, feature = "de_novo", pt.size.factor = 300)


###### signature expression in tumour regions

#LE
seurat_objects <- list(
  T243_RCTD_all_subpop_cnv = T243_RCTD_all_subpop_cnv,
  T241_RCTD_all_subpop_cnv = T241_RCTD_all_subpop_cnv,
  T242_RCTD_all_subpop_cnv = T242_RCTD_all_subpop_cnv,
  T256_RCTD_all_subpop_cnv = T256_RCTD_all_subpop_cnv,
  T260_RCTD_all_subpop_cnv = T260_RCTD_all_subpop_cnv,
  T262_RCTD_all_subpop_cnv = T262_RCTD_all_subpop_cnv,
  T266_RCTD_all_subpop_cnv = T266_RCTD_all_subpop_cnv,
  T269_RCTD_all_subpop_cnv = T269_RCTD_all_subpop_cnv,
  T270_RCTD_all_subpop_cnv = T270_RCTD_all_subpop_cnv,
  T275_RCTD_all_subpop_cnv = T275_RCTD_all_subpop_cnv,
  T296_RCTD_all_subpop_cnv = T296_RCTD_all_subpop_cnv,
  T304_RCTD_all_subpop_cnv = T304_RCTD_all_subpop_cnv,
  T334_RCTD_all_subpop_cnv = T334_RCTD_all_subpop_cnv,
  GBM_ZH881inf_RCTD_all_subpop_cnv = GBM_ZH881inf_RCTD_all_subpop_cnv,
  GBM_ZH881T1_RCTD_all_subpop_cnv = GBM_ZH881T1_RCTD_all_subpop_cnv,
  GBM_ZH916bulk_RCTD_all_subpop_cnv = GBM_ZH916bulk_RCTD_all_subpop_cnv,
  GBM_ZH916inf_RCTD_all_subpop_cnv = GBM_ZH916inf_RCTD_all_subpop_cnv,
  GBM_ZH1019inf_RCTD_all_subpop_cnv = GBM_ZH1019inf_RCTD_all_subpop_cnv,
  GBM_ZH8811Abulk_RCTD_all_subpop_cnv = GBM_ZH8811Abulk_RCTD_all_subpop_cnv,
  GBM_ZH8812bulk_RCTD_all_subpop_cnv = GBM_ZH8812bulk_RCTD_all_subpop_cnv
)

#CT
seurat_objects <- list(
  T243_RCTD_all_subpop_cnv = T243_RCTD_all_subpop_cnv,
  T242_RCTD_all_subpop_cnv = T242_RCTD_all_subpop_cnv,
  T248_RCTD_all_subpop_cnv = T248_RCTD_all_subpop_cnv,
  T251_RCTD_all_subpop_cnv = T251_RCTD_all_subpop_cnv,
  T255_RCTD_all_subpop_cnv = T255_RCTD_all_subpop_cnv,
  T256_RCTD_all_subpop_cnv = T256_RCTD_all_subpop_cnv,
  T259_RCTD_all_subpop_cnv = T259_RCTD_all_subpop_cnv,
  T260_RCTD_all_subpop_cnv = T260_RCTD_all_subpop_cnv,
  T262_RCTD_all_subpop_cnv = T262_RCTD_all_subpop_cnv,
  T266_RCTD_all_subpop_cnv = T266_RCTD_all_subpop_cnv,
  T268_RCTD_all_subpop_cnv = T268_RCTD_all_subpop_cnv,
  T269_RCTD_all_subpop_cnv = T269_RCTD_all_subpop_cnv,
  T275_RCTD_all_subpop_cnv = T275_RCTD_all_subpop_cnv,
  T304_RCTD_all_subpop_cnv = T304_RCTD_all_subpop_cnv,
  T313_RCTD_all_subpop_cnv = T313_RCTD_all_subpop_cnv,
  T334_RCTD_all_subpop_cnv = T334_RCTD_all_subpop_cnv,
  GBM_MGH258_RCTD_all_subpop_cnv = GBM_MGH258_RCTD_all_subpop_cnv,
  GBM_ZH881inf_RCTD_all_subpop_cnv = GBM_ZH881inf_RCTD_all_subpop_cnv,
  GBM_ZH881T1_RCTD_all_subpop_cnv = GBM_ZH881T1_RCTD_all_subpop_cnv,
  GBM_ZH916bulk_RCTD_all_subpop_cnv = GBM_ZH916bulk_RCTD_all_subpop_cnv,
  GBM_ZH916inf_RCTD_all_subpop_cnv = GBM_ZH916inf_RCTD_all_subpop_cnv,
  GBM_ZH916T1_RCTD_all_subpop_cnv = GBM_ZH916T1_RCTD_all_subpop_cnv,
  GBM_ZH1007inf_RCTD_all_subpop_cnv = GBM_ZH1007inf_RCTD_all_subpop_cnv,
  GBM_ZH1007nec_RCTD_all_subpop_cnv = GBM_ZH1007nec_RCTD_all_subpop_cnv,
  GBM_ZH1019inf_RCTD_all_subpop_cnv = GBM_ZH1019inf_RCTD_all_subpop_cnv,
  GBM_ZH1019T1_RCTD_all_subpop_cnv = GBM_ZH1019T1_RCTD_all_subpop_cnv,
  GBM_ZH8811Abulk_RCTD_all_subpop_cnv = GBM_ZH8811Abulk_RCTD_all_subpop_cnv,
  GBM_ZH8811Bbulk_RCTD_all_subpop_cnv = GBM_ZH8811Bbulk_RCTD_all_subpop_cnv,
  GBM_ZH8812bulk_RCTD_all_subpop_cnv = GBM_ZH8812bulk_RCTD_all_subpop_cnv
)

#CTmvp
seurat_objects <- list(
  T248_RCTD_all_subpop_cnv = T248_RCTD_all_subpop_cnv,
  T313_RCTD_all_subpop_cnv = T313_RCTD_all_subpop_cnv,
  GBM_ZH881T1_RCTD_all_subpop_cnv = GBM_ZH881T1_RCTD_all_subpop_cnv,
  GBM_ZH916bulk_RCTD_all_subpop_cnv = GBM_ZH916bulk_RCTD_all_subpop_cnv,
  GBM_ZH916T1_RCTD_all_subpop_cnv = GBM_ZH916T1_RCTD_all_subpop_cnv,
  GBM_ZH1007inf_RCTD_all_subpop_cnv = GBM_ZH1007inf_RCTD_all_subpop_cnv,
  GBM_ZH1007nec_RCTD_all_subpop_cnv = GBM_ZH1007nec_RCTD_all_subpop_cnv
)

#CTpan
seurat_objects <- list(
  T243_RCTD_all_subpop_cnv = T243_RCTD_all_subpop_cnv,
  T248_RCTD_all_subpop_cnv = T248_RCTD_all_subpop_cnv,
  T255_RCTD_all_subpop_cnv = T255_RCTD_all_subpop_cnv,
  T259_RCTD_all_subpop_cnv = T259_RCTD_all_subpop_cnv,
  T260_RCTD_all_subpop_cnv = T260_RCTD_all_subpop_cnv,
  T265_RCTD_all_subpop_cnv = T265_RCTD_all_subpop_cnv,
  T269_RCTD_all_subpop_cnv = T269_RCTD_all_subpop_cnv,
  T275_RCTD_all_subpop_cnv = T275_RCTD_all_subpop_cnv,
  T304_RCTD_all_subpop_cnv = T304_RCTD_all_subpop_cnv,
  T313_RCTD_all_subpop_cnv = T313_RCTD_all_subpop_cnv,
  T334_RCTD_all_subpop_cnv = T334_RCTD_all_subpop_cnv,
  GBM_ZH881T1_RCTD_all_subpop_cnv = GBM_ZH881T1_RCTD_all_subpop_cnv,
  GBM_ZH916bulk_RCTD_all_subpop_cnv = GBM_ZH916bulk_RCTD_all_subpop_cnv,
  GBM_ZH916T1_RCTD_all_subpop_cnv = GBM_ZH916T1_RCTD_all_subpop_cnv,
  GBM_ZH1007inf_RCTD_all_subpop_cnv = GBM_ZH1007inf_RCTD_all_subpop_cnv,
  GBM_ZH1007nec_RCTD_all_subpop_cnv = GBM_ZH1007nec_RCTD_all_subpop_cnv,
  GBM_ZH1019inf_RCTD_all_subpop_cnv = GBM_ZH1019inf_RCTD_all_subpop_cnv,
  GBM_ZH1019T1_RCTD_all_subpop_cnv = GBM_ZH1019T1_RCTD_all_subpop_cnv,
  GBM_ZH8811Abulk_RCTD_all_subpop_cnv = GBM_ZH8811Abulk_RCTD_all_subpop_cnv,
  GBM_ZH8811Bbulk_RCTD_all_subpop_cnv = GBM_ZH8811Bbulk_RCTD_all_subpop_cnv
)


# Define the metadata column of interest
metadata_column <- "de_novo"  # Replace with the actual metadata column

# Define the region of interest in Tumour_region_SpotIDs
region_of_interest <- "CTmvp_genes"  #(make sure to select corresponding Seurat object list)

# Initialize an empty dataframe to store results
results <- data.frame(
  Sample = character(),
  Region = character(),
  AvgScore_Region = numeric(),
  AvgScore_non_Region = numeric(),
  log2FC = numeric(),
  stringsAsFactors = FALSE
)

# Set a seed for reproducibility
set.seed(1234)

# Iterate over each Seurat object
for (sample_name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[sample_name]]
  
  # Ensure the metadata column exists
  if (!(metadata_column %in% colnames(seurat_obj@meta.data))) {
    message(paste("Skipping", sample_name, "- Metadata column not found"))
    next
  }
  
  # Get all spots labeled as the chosen region of interest
  region_spots <- WhichCells(seurat_obj, expression = Tumour_region_SpotIDs == region_of_interest)
  
  # Get all non-selected spots
  non_region_spots <- setdiff(Cells(seurat_obj), region_spots)
  
  # Check if there are any region-specific spots
  if (length(region_spots) == 0 || length(non_region_spots) == 0) {
    message(paste("Skipping", sample_name, "- No spots found for", region_of_interest))
    next
  }
  
  # Determine the number of spots to use (minimum of region or non-region)
  num_spots <- min(length(region_spots), length(non_region_spots))
  
  # Set seed before sampling to ensure reproducibility
  set.seed(1234)
  selected_region_spots <- sample(region_spots, num_spots)
  
  set.seed(1234)
  selected_non_region_spots <- sample(non_region_spots, num_spots)
  
  # Extract metadata score for selected spots
  score_region <- FetchData(seurat_obj, vars = metadata_column, cells = selected_region_spots)
  avg_score_region <- mean(score_region[[metadata_column]], na.rm = TRUE)
  
  # Extract metadata score for randomly selected non-region spots
  score_non_region <- FetchData(seurat_obj, vars = metadata_column, cells = selected_non_region_spots)
  avg_score_non_region <- mean(score_non_region[[metadata_column]], na.rm = TRUE)
  
  # Calculate log2 Fold Change (log2FC)
  log2_FC <- log2((avg_score_region + 1) / (avg_score_non_region + 1))
  
  # Store results
  results <- rbind(
    results, 
    data.frame(
      Sample = sample_name, 
      Region = region_of_interest,
      AvgScore_Region = avg_score_region, 
      AvgScore_non_Region = avg_score_non_region,
      log2FC = log2_FC
    )
  )
}

# Display results
print(results)

###### ANOVA for spatial enrichment

# Required libraries
library(dplyr)
library(stats)

# Load your data
data <- read.csv("file location for spatial enrichment calulations csv")

# Convert 'Region' to a factor
data$Region <- as.factor(data$Region)

# Get the list of genes (excluding 'Sample' and 'Region' columns)
gene_list <- colnames(data)[3:ncol(data)]

# Initialize a list to store results
anova_results <- list()
tukey_results <- list()
adjusted_pvalues <- list()

# Perform ANOVA and post-hoc tests
for (gene in gene_list) {
  formula <- as.formula(paste(gene, "~ Region"))
  model <- aov(formula, data = data)
  
  # Store ANOVA result
  anova_summary <- summary(model)
  p_value <- anova_summary[[1]]$`Pr(>F)`[1]  # Extract p-value from ANOVA
  anova_results[[gene]] <- p_value
  
  # Post-hoc Tukey test
  tukey_test <- TukeyHSD(model)
  tukey_results[[gene]] <- tukey_test$Region
  
  # Store raw p-values from post-hoc tests
  raw_pvalues <- tukey_test$Region[, 4]  # Extract post-hoc p-values
  
  # Apply Benjamini-Hochberg correction
  adjusted_pvalues[[gene]] <- p.adjust(raw_pvalues, method = "BH")
}

# Convert results to data frames for easier visualization
anova_df <- data.frame(Gene = names(anova_results), ANOVA_pvalue = unlist(anova_results))
tukey_df <- do.call(rbind, lapply(names(tukey_results), function(gene) {
  data.frame(Gene = gene, Comparison = rownames(tukey_results[[gene]]), 
             Tukey_pvalue = tukey_results[[gene]][, 4],
             BH_adjusted_pvalue = adjusted_pvalues[[gene]])
}))

# Print results
print("ANOVA p-values:")
print(anova_df)

print("Tukey post-hoc results with BH correction:")
print(tukey_df)


######## Generate Spatial plots

# Plot the Genesofinterest score (optional)
SpatialFeaturePlot(GBM_ZH8811Bbulk_RCTD_all_subpop_cnv, feature = "WDR77", pt.size.factor = 3.5)
SpatialDimPlot(GBM_ZH8811Bbulk_RCTD_all_subpop_cnv, group.by = "Tumour_region_SpotIDs", pt.size.factor = 3.5)

