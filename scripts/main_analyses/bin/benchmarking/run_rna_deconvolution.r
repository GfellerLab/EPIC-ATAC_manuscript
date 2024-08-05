

# Libraries ---------------------------------------------------------------
library(EPIC)
library(dplyr)
library(tidyr)


# Functions ---------------------------------------------------------------

load_bulk_proportion = function(bulks_data_rna, dataset_name, renorm = F){
  bulk_proportions = bulks_data_rna$obs
  bulk_proportions = bulk_proportions[, c(1:(ncol(bulk_proportions) - 1)), drop = F]
  head(bulk_proportions)
  bulk_proportions = as.data.frame(t(bulk_proportions))
  
  if(dataset_name %in% c("PBMC_multiome_simulated")){
    bulk_proportions$CD4_Tcells = rowSums(bulk_proportions[, c("Memory_CD4", "Naive_CD4", "Treg")])
    bulk_proportions = bulk_proportions[, -which(colnames(bulk_proportions) %in% c("Memory_CD4", "Naive_CD4", "Treg"))]
    bulk_proportions$CD8_Tcells = rowSums(bulk_proportions[, c("Naive_CD8_Tcells", "Non_Naive_CD8_Tcells")])
    bulk_proportions = bulk_proportions[, -which(colnames(bulk_proportions) %in% c("Naive_CD8_Tcells", "Non_Naive_CD8_Tcells"))]
  } 
  
  if(renorm){
    bulk_proportions = bulk_proportions[,-which(colnames(bulk_proportions) %in% c("Endothelial", "Fibroblasts", "Other"))] 
    bulk_proportions = bulk_proportions/rowSums(bulk_proportions)
  }  
  return(bulk_proportions)
} 

# Parameters --------------------------------------------------------------
bulk_path = as.character(commandArgs(TRUE)[1])

dataset_name =  gsub("_bulk_RNA.rdata", "", basename(bulk_path))
load(bulk_path)

bulk_proportions = load_bulk_proportion(bulks_data_rna, dataset_name)

# Run EPIC ----------------------------------------------
epic_rna = EPIC(
  bulk = bulks_data_rna$bulk,
  reference = ifelse(grepl("PBMC", dataset_name), "BRef", "TRef"),
  withOtherCells = ifelse(grepl("PBMC", dataset_name), F, T)
)

epic_rna = as.data.frame(epic_rna$cellFractions)
head(epic_rna)

colnames(epic_rna) <- plyr::revalue(colnames(epic_rna), c("CAFs" = "Fibroblasts", "otherCells" = "Other")) #"NKcells" = "NK",
if(!grepl("PBMC", dataset_name)){
  epic_rna$Tcells <- epic_rna$CD4_Tcells + epic_rna$CD8_Tcells
}
epic_rna = epic_rna[, which(colnames(epic_rna) %in% colnames(bulk_proportions))]
bulk_proportions = bulk_proportions[, which(colnames(bulk_proportions) %in%  colnames(epic_rna))] 

summary_deconv = data.frame(
  cell_type = as.vector(sapply(1:ncol(bulk_proportions), function(i) rep(colnames(bulk_proportions)[i], nrow(bulk_proportions)))),
  sample = rep(rownames(epic_rna), ncol(epic_rna)),
  estimate = as.vector(sapply(colnames(bulk_proportions), function(i) epic_rna[, i])),
  method = "EPIC_RNA",
  true_fraction = as.vector(sapply(1:ncol(bulk_proportions), function(i) bulk_proportions[, i])),
  bulk_name = rep(dataset_name, ncol(epic_rna))
)  


write.table(summary_deconv, file = paste0(dataset_name, "_RNA_deconvolution_summary.txt"), quote = F, row.names = F, col.names = T, sep = "\t")



