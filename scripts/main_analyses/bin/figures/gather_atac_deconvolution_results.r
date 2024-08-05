# Libraries ---------------------------------------------------------------
library(readxl)
library(ggplot2)
library(Metrics)

# Functions ---------------------------------------------------------------

get_pred_truth = function(deconv_res, bulk_proportions, cells_matching, method = "EPIC"){
  head(deconv_res)
  
  deconv_comp = data.frame(sample = vector(),
                           obs = vector(),
                           cell = vector(),
                           pred = vector())
  
  for(cell_type in unique(cells_matching$Grouping_name_bulk)){
    
    matching_rows = cells_matching[which(cells_matching[, "Grouping_name_bulk"] == cell_type), , drop = F] 
    
    prediction_columns = unique(matching_rows[, paste0("Cell_types_pred_",method)])
    prediction_columns = unique(as.vector(unlist(sapply(prediction_columns, function(i) unlist(strsplit(i,","))))))
     
    if(!"NA" %in% prediction_columns ){
      if(sum(prediction_columns %in% colnames(deconv_res)) != 0 ){
        cell_type_bulk_proportions = rowSums(bulk_proportions[, matching_rows[, "Cell_type_bulk"], drop = F])
        cell_type_predictions = rowSums(deconv_res[, prediction_columns , drop =F])
      } else{
        cell_type_bulk_proportions = rowSums(bulk_proportions[, matching_rows[, "Cell_type_bulk"], drop = F])
        cell_type_predictions = rep(NA, nrow(deconv_res))
      } 
      
    } else{
      cell_type_bulk_proportions = rowSums(bulk_proportions[, matching_rows[, "Cell_type_bulk"], drop = F])
      cell_type_predictions = rep(0, nrow(deconv_res))
    } 
    
    deconv_comp = rbind(
      deconv_comp,
      data.frame(
        sample = rownames(bulk_proportions),
        obs = cell_type_bulk_proportions,
        cell = rep(cell_type, nrow(bulk_proportions)),
        pred = cell_type_predictions
      )
    )
    
  } 
  
  considered_cell_types = unique(cells_matching[, paste0("Cell_types_pred_",method)])
  considered_cell_types = unique(as.vector(unlist(sapply(considered_cell_types, function(i) unlist(strsplit(i,","))))))

  cell_types_absent_in_bulk = colnames(deconv_res)[!colnames(deconv_res) %in% considered_cell_types]
  if (length(cell_types_absent_in_bulk) != 0) {
    for (cell_type in cell_types_absent_in_bulk) {
      deconv_comp = rbind(
        deconv_comp,
        data.frame(
          sample = rownames(bulk_proportions),
          obs = rep(0, nrow(deconv_res)),
          cell = rep(cell_type, nrow(bulk_proportions)),
          pred = deconv_res[, cell_type]
        )
      )
    }
  } 
  
  return(deconv_comp)
}


# Parameters --------------------------------------------------------------
benchmarking_output_path = as.character(commandArgs(TRUE)[1])
path_cell_matching = as.character(commandArgs(TRUE)[2])
with_subtypes = as.logical(commandArgs(TRUE)[3])

if(with_subtypes){
  bulk_names = c("BCC_Satpathy", "PBMC_Granja", "PBMC_Satpathy", "PBMC_experiment", "PBMC_multiome") 
}else{
  bulk_names = c("BCC_Satpathy", "Gynecological_cancers", "PBMC_Granja", "PBMC_Satpathy", "PBMC_experiment", "PBMC_multiome", "PBMC_multiome_simulated","HNSCC","BRCA","CRC","CESC","SKCM","PDAC", "OV") 
}
summary_predictions = data.frame(cell_type = vector(),
                                 sample = vector(),
                                 estimate = vector(),
                                 method = vector(),
                                 true_fraction = vector(),
                                 bulk_name = vector())      
for(bulk in bulk_names){
  
  # Load true proportions ---------------------------------------------------
  bulk_proportions = read.table(paste0(benchmarking_output_path, "intermediate_files/", bulk, "/", bulk, "_bulk_proportions.txt"), header = T)
  rownames(bulk_proportions) = bulk_proportions$cell_type
  bulk_proportions = bulk_proportions[, c(1:(ncol(bulk_proportions) - 2)), drop = F]
  head(bulk_proportions)
  bulk_proportions = as.data.frame(t(bulk_proportions))
  
  # Read cells matching -----------------------------------------------------
  cells_matching = as.data.frame(read_excel(paste0(path_cell_matching)))
  head(cells_matching)
  cells_matching = cells_matching[which(cells_matching$Data_type == paste0("Bulk_", bulk)), ]
  
  # EPIC results ------------------------------------------------------------
  deconvolution_output_path = list.files(paste0(benchmarking_output_path, "R_methods/"), paste0(bulk,"_EPIC_predictions.txt"), recursive = T, full.names = T)
   
  deconv_res = read.table(deconvolution_output_path, header = T)
  deconv_res = deconv_res[, c(1:(ncol(deconv_res) - 1)), drop = F]
  
  deconv_comp = get_pred_truth(deconv_res, bulk_proportions, cells_matching, method = "withOtherCells")
  summary_predictions = rbind(summary_predictions, data.frame(cell_type = deconv_comp$cell,
                                                              sample = deconv_comp$sample,
                                                              estimate = deconv_comp$pred,
                                                              method = "EPIC-ATAC",
                                                              true_fraction = deconv_comp$obs,
                                                              bulk_name = rep(bulk, nrow(deconv_comp))))   
  
  # Quantiseq results ------------------------------------------------------------
  deconvolution_output_path = list.files(paste0(benchmarking_output_path, "R_methods/"), paste0(bulk,"_quantiseq_predictions"), recursive = T, full.names = T)
  deconv_res = read.table(deconvolution_output_path, header = T)
  deconv_res = deconv_res[, c(1:(ncol(deconv_res) - 1)), drop = F]
  colnames(deconv_res)[ncol(deconv_res)] = "otherCells"  
  
  deconv_comp = get_pred_truth(deconv_res, bulk_proportions, cells_matching, method = "withOtherCells")
  summary_predictions = rbind(summary_predictions, data.frame(cell_type = deconv_comp$cell,
                                                              sample = deconv_comp$sample,
                                                              estimate = deconv_comp$pred,
                                                              method = "quanTIseq",
                                                              true_fraction = deconv_comp$obs,
                                                              bulk_name = rep(bulk, nrow(deconv_comp))))   
  
  # ABIS results ------------------------------------------------------------
  deconvolution_output_path = list.files(paste0(benchmarking_output_path, "R_methods/"), paste0(bulk,"_ABIS_predictions"), recursive = T, full.names = T)
  deconv_res = read.table(deconvolution_output_path, header = T)
  deconv_res = deconv_res[, c(1:(ncol(deconv_res) - 1)), drop = F]

  deconv_comp = get_pred_truth(deconv_res, bulk_proportions, cells_matching, method = "noOtherCells")
  summary_predictions = rbind(summary_predictions, data.frame(cell_type = deconv_comp$cell,
                                                              sample = deconv_comp$sample,
                                                              estimate = deconv_comp$pred,
                                                              method = "ABIS",
                                                              true_fraction = deconv_comp$obs,
                                                              bulk_name = rep(bulk, nrow(deconv_comp))))   
  
  # MCPcounter results ------------------------------------------------------------
  deconvolution_output_path = list.files(paste0(benchmarking_output_path, "R_methods/"), paste0(bulk,"_MCPcounter_predictions"), recursive = T, full.names = T)
  deconv_res = read.table(deconvolution_output_path, header = T)
  deconv_res = deconv_res[, c(1:(ncol(deconv_res) - 1)), drop = F]
  
  if(grepl("PBMC", bulk)){
    deconv_res = deconv_res[ ,which(!colnames(deconv_res) %in% c("Fibroblasts", "Endothelial", "Macrophages"))] 
  } else{
    deconv_res = deconv_res[ ,which(colnames(deconv_res) != "Monocytes")] 
  } 
  deconv_comp = get_pred_truth(deconv_res, bulk_proportions, cells_matching, method = "noOtherCells")
  summary_predictions = rbind(summary_predictions, data.frame(cell_type = deconv_comp$cell,
                                                              sample = deconv_comp$sample,
                                                              estimate = deconv_comp$pred,
                                                              method = "MCPcounter",
                                                              true_fraction = deconv_comp$obs,
                                                              bulk_name = rep(bulk, nrow(deconv_comp)))) 
  # Deconpeaker results with original markers ------------------------------------------------------------
  deconvolution_output_path = list.files(paste0(benchmarking_output_path, "DeconPeaker_original_markers/deconvolution/"), paste0(bulk,"_deconPeaker-Results_original_markers.xls"), recursive = T, full.names = T)
  deconv_res = read.table(deconvolution_output_path, header = T)
  deconv_res = deconv_res[, c(1:(ncol(deconv_res) - 3)), drop = F]
  
  deconv_comp = get_pred_truth(deconv_res, bulk_proportions, cells_matching, method = "deconpeaker_OM")
  summary_predictions = rbind(summary_predictions, data.frame(cell_type = deconv_comp$cell,
                                                              sample = deconv_comp$sample,
                                                              estimate = deconv_comp$pred,
                                                              method = "DeconPeaker_originalMarkers",
                                                              true_fraction = deconv_comp$obs,
                                                              bulk_name = rep(bulk, nrow(deconv_comp)))) 
  
 # Deconpeaker results with our markers ------------------------------------------------------------
  deconvolution_output_path = list.files(paste0(benchmarking_output_path, "DeconPeaker/deconvolution/"), paste0(bulk,"_deconPeaker-Results_epic_markers.xls"), recursive = T, full.names = T)
  deconv_res = read.table(deconvolution_output_path, header = T)
  deconv_res = deconv_res[, c(1:(ncol(deconv_res) - 3)), drop = F]
  
  deconv_comp = get_pred_truth(deconv_res, bulk_proportions, cells_matching, method = "noOtherCells")
  summary_predictions = rbind(summary_predictions, data.frame(cell_type = deconv_comp$cell,
                                                              sample = deconv_comp$sample,
                                                              estimate = deconv_comp$pred,
                                                              method = "DeconPeaker_ourMarkers",
                                                              true_fraction = deconv_comp$obs,
                                                              bulk_name = rep(bulk, nrow(deconv_comp))))    
  
  # Deconpeaker results with the reference generated by Deconpeaker ------------------------------------------------------------
  deconvolution_output_path = list.files(paste0(benchmarking_output_path, "DeconPeaker_newRef/deconvolution/"), paste0(bulk,"_deconPeaker-Results_newRef.xls"), recursive = T, full.names = T)
  deconv_res=read.table(deconvolution_output_path,header = T)
  deconv_res=deconv_res[,c(1:(ncol(deconv_res)-3)),drop=F]
  
  deconv_comp = get_pred_truth(deconv_res, bulk_proportions, cells_matching, method = "noOtherCells")
  summary_predictions = rbind(summary_predictions, data.frame(cell_type = deconv_comp$cell,
                                                              sample = deconv_comp$sample,
                                                              estimate = deconv_comp$pred,
                                                              method = "DeconPeaker-Custom",
                                                              true_fraction = deconv_comp$obs,
                                                              bulk_name = rep(bulk, nrow(deconv_comp)))) 
  
  # CIBERSORTx results with our markers ------------------------------------------------------------
  deconvolution_output_path = list.files(paste0(benchmarking_output_path, "CIBERSORTx_ourMarkers"), paste0(bulk,"_CIBERSORTx_noAbs_Results_ourMarkers.txt"), recursive = T, full.names = T)
  deconv_res=read.table(deconvolution_output_path,header = T, sep = "\t")
  deconv_res=deconv_res[,c(2:(ncol(deconv_res)-3)),drop=F]

  deconv_comp = get_pred_truth(deconv_res, bulk_proportions, cells_matching, method = "noOtherCells")
  summary_predictions = rbind(summary_predictions, data.frame(cell_type = deconv_comp$cell,
                                                              sample = deconv_comp$sample,
                                                              estimate = deconv_comp$pred,
                                                              method = "CIBERSORTx_ourMarkers",
                                                              true_fraction = deconv_comp$obs,
                                                              bulk_name = rep(bulk, nrow(deconv_comp))))

  # CIBERSORTx results with the reference generated by CIBERSORTx ------------------------------------------------------------
  deconvolution_output_path = list.files(paste0(benchmarking_output_path, "CIBERSORTx_newRef"), paste0(bulk,"_CIBERSORTx_noAbs_Results_newRef.txt"), recursive = T, full.names = T)
  deconv_res = read.table(deconvolution_output_path, header = T, sep = "\t")
  deconv_res = deconv_res[, c(2:(ncol(deconv_res) - 3)), drop = F]

  deconv_comp = get_pred_truth(deconv_res, bulk_proportions, cells_matching, method = "noOtherCells")
  summary_predictions = rbind(summary_predictions, data.frame(cell_type = deconv_comp$cell,
                                                              sample = deconv_comp$sample,
                                                              estimate = deconv_comp$pred,
                                                              method = "CIBERSORTx-Custom",
                                                              true_fraction = deconv_comp$obs,
                                                              bulk_name = rep(bulk, nrow(deconv_comp))))
} 



# Renaming before figure generations 
summary_predictions$bulk_name_original = summary_predictions$bulk_name

summary_predictions$bulk_name = plyr::revalue(summary_predictions$bulk_name, c("BCC_Satpathy"="Basal cell carcinoma"))
summary_predictions$sample[which(summary_predictions$bulk_name == "PBMC_Granja")] ="Granja"
summary_predictions$sample[which(summary_predictions$bulk_name == "PBMC_Satpathy")] ="Satpathy"
summary_predictions$sample[which(summary_predictions$bulk_name == "PBMC_multiome")] ="10X_multiome"
summary_predictions$bulk_name[which(summary_predictions$bulk_name %in% c("HNSCC","BRCA","CRC","CESC","SKCM","PDAC", "OV", "GBM", "UCEC", "CEAD"))] = "HTAN"

summary_predictions$bulk_name[which(summary_predictions$bulk_name %in% c("PBMC_Granja","PBMC_Satpathy","PBMC_multiome"))]="PBMC pseudobulk"
summary_predictions$bulk_name <- plyr::revalue(summary_predictions$bulk_name, c("PBMC_experiment" = "PBMC experiment",
                                                                                "Gynecological_cancers" = "Gynecological cancers"))
summary_predictions$cell_type <- plyr::revalue(summary_predictions$cell_type, c("Other" = "Uncharacterized", 
                                                                                "Myeloid" = "DCs + Macrophages", 
                                                                                "Non_Naive_CD4_Tcells" = "Memory + Helper CD4",
                                                                                "Naive_CD4_Tcells" = "Naive CD4", 
                                                                                "Non_Naive_CD8_Tcells" = "Non Naive CD8", 
                                                                                "Naive_CD8_Tcells" = "Naive CD8",
                                                                                "NK" = "NKcells",
                                                                                "DCs" = "Dendritic",
                                                                                "NK_T_cells" = "NK + T cells"))

# Save the predictions summary --------------------------------
write.table(summary_predictions, file = paste0("summary_atac_predictions.txt"), quote = F, row.names = F, col.names = T, sep = "\t")


# Renormalize deconpeaker and cibersortx true proportions when otherCells = T ----------------------
for(bulk in c("Basal cell carcinoma","Gynecological cancers", "HTAN")){ 
  
  if(bulk %in% summary_predictions$bulk_name){
    for(method in unique(summary_predictions$method)){
      summary_predictions = summary_predictions[-which(summary_predictions$method == method &
                                                         summary_predictions$bulk_name == bulk &
                                                         summary_predictions$cell_type == "Uncharacterized"), ] 
      
      samples = unique(summary_predictions$sample[which(summary_predictions$method == method & summary_predictions$bulk_name == bulk)] )
      for(sample in samples){
        prop1 = summary_predictions$estimate[which(summary_predictions$method == method & 
                                                     summary_predictions$bulk_name == bulk & 
                                                     summary_predictions$sample == sample )] 
        prop2 = (prop1)/sum(prop1,na.rm = T)
        summary_predictions$estimate[which(summary_predictions$method == method & summary_predictions$bulk_name == bulk & summary_predictions$sample == sample)] = prop2
        prop1 = summary_predictions$true_fraction[which(summary_predictions$method == method & 
                                                          summary_predictions$bulk_name == bulk & 
                                                          summary_predictions$sample == sample )] 
        prop2 = (prop1)/sum(prop1, na.rm = T)
        summary_predictions$true_fraction[which(summary_predictions$method == method & summary_predictions$bulk_name == bulk & summary_predictions$sample == sample)] = prop2
      }  
    } 
  }
  
} 


# Save the predictions summary --------------------------------
write.table(summary_predictions, file = paste0("summary_atac_predictions-renorm.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

