############## Heatmap for KDM gene family scRNAseq data ###########

#Load necessary libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)

###load datasets

#Ebert
obj_Ebert_Mal <- readRDS("file path/Ebert.rds")
Idents(obj_Ebert) <- "cell_type"
obj_Ebert_Mal <- subset(obj_Ebert, idents = c("malignant"))

#LeBlanc
obj_LeBlanc <- readRDS("file path/LeBlanc.rds")
Idents(obj_LeBlanc) <- "sample"
table(obj_LeBlanc$sample)
obj_LeBlanc <- subset(obj_LeBlanc, idents = c("JK124_reg1_tis_1", "JK124_reg1_tis_2", "JK124_reg2_tis_1", "JK124_reg2_tis_2", "JK125_reg1_tis_1.1", "JK125_reg1_tis_1.2", "JK125_reg2_tis_1", "JK125_reg2_tis_2_r1", "JK126_reg1_tis_1.1", "JK126_reg1_tis_1.2", "JK126_reg2_tis_1", "JK134_reg1_tis_1", "JK134_reg2_tis_1", "JK136_reg1_tis_1", "JK136_reg2_tis_1", "JK136_reg2_tis_2_br", "JK142_reg1_tis_1", "JK142_reg2_tis_1", "JK142_reg2_tis_2.1_br", "JK142_reg2_tis_2.2_br", "JK152_reg1_tis_1", "JK152_reg2_tis_1", "JK153_reg1_tis_1", "JK153_reg2_tis_1", "JK156_reg1_tis_1", "JK156_reg2_tis_1", "JK156_reg2_tis_2_br", "JK163_reg1_tis_1", "JK163_reg2_tis_1")) #subset primary GBM tissues
Idents(obj_LeBlanc) <- "cell_type"
obj_LeBlanc_Mal <- subset(obj_LeBlanc, idents = c("malignant"))
table(obj_LeBlanc_Mal$patient)

#Neftel
obj_Neftel <- readRDS("file path/Neftel.rds")
Idents(obj_Neftel) <- "GBMType"
obj_Neftel <- subset(obj_Neftel, idents = c("Adult"))
Idents(obj_Neftel) <- "CellAssignment"
obj_Neftel_Mal <- subset(obj_Neftel, idents = c("Malignant"))
table(obj_Neftel_Mal$Sample)


#Abdelfattah
obj_Abdelfattah <- readRDS("file path/Abdelfattah.rds")
Idents(obj_Abdelfattah) <- "Type"
obj_Abdelfattah <- subset(obj_Abdelfattah, idents = c("GBM"))
Idents(obj_Abdelfattah) <- "Assignment"
obj_Abdelfattah_Mal <- subset(obj_Abdelfattah, idents = c("Glioma"))
table(obj_Abdelfattah_Mal$Patient)
Idents(obj_Abdelfattah_Mal) <- "Patient"
obj_Abdelfattah_Mal <- subset(obj_Abdelfattah_Mal, idents = c("ndGBM-01", "ndGBM-02", "ndGBM-04", "ndGBM-05", "ndGBM-06", "ndGBM-07", "ndGBM-08", "ndGBM-09")) #subset primary GBM tissues


#Chen
obj_Chen <- readRDS("file path/Chen.rds")
Idents(obj_Chen) <- "orig.ident"
obj_Chen <- subset(obj_Chen, idents = c("PDC001", "PJ052", "PJ053", "PW016-703", "PW017-703", "PW032-710", "PW035-710", "PW039-705"))  #subset primary GBM tissues
Idents(obj_Chen) <- "cellkb"
obj_Chen_Mal <- subset(obj_Chen, idents = c("Tumor"))
table(obj_Chen_Mal$orig.ident)


# Couturier ####
obj_Couturier <- readRDS("file path/Couturier.rds")
Idents(obj_Couturier) <- "orig.ident"
obj_Couturier <- subset(obj_Couturier, idents = c("BT333", "BT338_1of2", "BT338_2of2", "BT346", "BT363_1of2", "BT363_2of2", "BT364_1of2", "BT364_2of2", "BT368", "BT389", "BT390", "BT397_1of2", "BT397_2of2", "BT400", "BT402", "BT407", "BT409"))  #subset primary GBM tissues
Idents(obj_Couturier) <- "cellkb"
obj_Couturier_Mal <- subset(obj_Couturier, idents = c("Tumor"))
table(obj_Couturier_Mal$orig.ident)
obj_Couturier_Mal$Sample <- gsub("(_1of2|_2of2)$", "", obj_Couturier_Mal$orig.ident) #make a new sample ID metadata column to merge samples split in two
table(obj_Couturier_Mal$Sample)

######################### Calculate z scores within each Patient ############################

library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)

# Genes of interest (KDM family)
genes_of_interest <- c("KDM1A", "KDM1B", "KDM2A", "KDM2B", "KDM3A", "KDM3B", "JMJD1C", "KDM4A", 
                       "KDM4B", "KDM4C", "KDM4D", "KDM4E", "KDM4F", "KDM5A", "KDM5B", "KDM5C", "KDM5D", "KDM6A", 
                       "KDM6B", "KDM7A", "PHF8", "PHF2", "KDM8")

# Seurat objects
seurat_list <- list(
  Ebert = obj_Ebert_Mal,
  LeBlanc = obj_LeBlanc_Mal,
  Neftel = obj_Neftel_Mal,
  Abdelfattah = obj_Abdelfattah_Mal,
  Chen = obj_Chen_Mal,
  Couturier = obj_Couturier_Mal
)

# Corresponding metadata column names for patient IDs
patient_id_cols <- list(
  Ebert = "orig.ident",
  LeBlanc = "patient",
  Neftel = "Sample",
  Abdelfattah = "Patient",
  Chen = "orig.ident",
  Couturier = "Sample"
)

# Compute z-scores "within" each patient (i.e. patient wise z score)
zscore_by_patient <- function(seurat_obj, dataset_name) {
  id_col <- patient_id_cols[[dataset_name]]
  present_genes <- intersect(genes_of_interest, rownames(seurat_obj))
  
  if (length(present_genes) == 0) {
    warning(paste0("No genes of interest found in ", dataset_name))
    return(NULL)
  }
  
  expr_data <- GetAssayData(seurat_obj, layer = "data")[present_genes, , drop = FALSE]
  cell_metadata <- seurat_obj@meta.data
  patient_ids <- cell_metadata[[id_col]]
  names(patient_ids) <- colnames(seurat_obj)
  
  avg_expr_by_patient <- as.data.frame(expr_data) %>%
    t() %>%
    as.data.frame() %>%
    mutate(patient_id = patient_ids) %>%
    group_by(patient_id) %>%
    summarise(across(all_of(present_genes), \(x) mean(x, na.rm = TRUE))) %>%
    pivot_longer(-patient_id, names_to = "gene", values_to = "avg_expression") %>%
    group_by(patient_id) %>%
    mutate(z_score = scale(avg_expression)[, 1]) %>%
    ungroup() %>%
    mutate(dataset = dataset_name)
  
  return(avg_expr_by_patient)
}

# Apply function to all datasets
zscore_all <- map2_dfr(seurat_list, names(seurat_list), zscore_by_patient)

# transform to wide matrix: genes as rows, dataset_patient as columns
zscore_matrix <- zscore_all %>%
  select(gene, patient_id, dataset, z_score) %>%
  unite("sample_id", dataset, patient_id, sep = "_") %>%
  pivot_wider(names_from = sample_id, values_from = z_score) %>%
  arrange(gene)

# save file for import into https://software.broadinstitute.org/morpheus/
write.csv(zscore_matrix, "file path/zscore_matrix_within_patient_KDM_family.csv")

